#!/usr/bin/env python3
"""Compute TBLite forces for extxyz crystal cells and write force-annotated extxyz files."""

import argparse
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path


def _set_single_thread_env() -> None:
    """Keep each worker to one native thread; parallelism is across structures."""
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[variable] = "1"


_set_single_thread_env()

import pandas as pd
import numpy as np
from ase.constraints import FixAtoms
from ase.io import read, write
from ase.optimize import FIRE
from tblite.ase import TBLite


DEFAULT_WORKERS = os.cpu_count() or 1
SUMMARY_PATH = "tblite_force_summary.csv"
DEFAULT_MAX_ATOMS = 300
DEFAULT_H_FMAX = 0.05
DEFAULT_H_STEPS = 100


def _clear_previous_outputs(output_dir: Path) -> None:
    for path in output_dir.glob("*.xyz"):
        path.unlink()
    summary_path = output_dir / SUMMARY_PATH
    if summary_path.exists():
        summary_path.unlink()


def _tblite_calculator() -> TBLite:
    return TBLite(
        method="GFN2-xTB",
        charge=0,
        multiplicity=1,
        verbosity=0,
    )


def _force_summary_row(
    input_path: Path,
    out_path: Path,
    status: str,
    n_atoms: int | None,
    max_force_ev_per_ang: float | None,
    error: str = "",
) -> dict:
    return {
        "mol_id": input_path.stem,
        "input_path": str(input_path),
        "extxyz_path": str(out_path),
        "status": status,
        "n_atoms": n_atoms,
        "max_force_ev_per_ang": max_force_ev_per_ang,
        "error": error,
    }


def _write_force_extxyz(atoms, out_path: Path) -> None:
    write(
        out_path,
        atoms,
        format="extxyz",
        columns=["symbols", "positions", "forces"],
        write_results=True,
    )


def _calculate_and_write_forces(atoms, input_path: Path, out_path: Path) -> dict:
    atoms.calc = _tblite_calculator()
    forces = atoms.get_forces()
    _write_force_extxyz(atoms, out_path)
    return _force_summary_row(
        input_path=input_path,
        out_path=out_path,
        status="ok",
        n_atoms=len(atoms),
        max_force_ev_per_ang=float(np.linalg.norm(forces, axis=1).max()),
    )


def _optimized_summary_row(
    base_row: dict,
    h_atoms: int | None,
    optimization_converged: bool | None,
    optimization_steps: int | None,
    max_h_force_ev_per_ang: float | None,
) -> dict:
    return {
        **base_row,
        "h_atoms": h_atoms,
        "optimization_converged": optimization_converged,
        "optimization_steps": optimization_steps,
        "max_h_force_ev_per_ang": max_h_force_ev_per_ang,
    }


def _optimize_hydrogens_and_write_forces(
    atoms,
    input_path: Path,
    out_path: Path,
    h_fmax: float,
    h_steps: int,
) -> dict:
    h_mask = np.asarray(atoms.symbols == "H")
    h_atoms = int(h_mask.sum())
    atoms.calc = _tblite_calculator()

    if h_atoms:
        atoms.set_constraint(FixAtoms(mask=~h_mask))
        optimizer = FIRE(atoms, logfile=None)
        optimization_converged = bool(optimizer.run(fmax=h_fmax, steps=h_steps))
        optimization_steps = optimizer.nsteps
    else:
        optimization_converged = True
        optimization_steps = 0

    # ASE optimizers use Atoms.get_forces() with constraints applied. For
    # reporting, request raw calculator forces so fixed heavy atoms are not
    # written with zero force vectors.
    forces = atoms.get_forces(apply_constraint=False)
    _write_force_extxyz(atoms, out_path)

    force_norms = np.linalg.norm(forces, axis=1)
    max_h_force = float(force_norms[h_mask].max()) if h_atoms else None
    return _optimized_summary_row(
        _force_summary_row(
            input_path=input_path,
            out_path=out_path,
            status="ok",
            n_atoms=len(atoms),
            max_force_ev_per_ang=float(force_norms.max()),
        ),
        h_atoms=h_atoms,
        optimization_converged=optimization_converged,
        optimization_steps=optimization_steps,
        max_h_force_ev_per_ang=max_h_force,
    )


def _failed_row(input_path: Path, out_path: Path, n_atoms: int | None, exc: Exception) -> dict:
    return _force_summary_row(
        input_path=input_path,
        out_path=out_path,
        status="failed",
        n_atoms=n_atoms,
        max_force_ev_per_ang=None,
        error=f"{type(exc).__name__}: {exc}",
    )


def _compute_one(
    input_path: Path,
    output_dir: Path,
    optimized_output_dir: Path,
    max_atoms: int | None,
    h_fmax: float,
    h_steps: int,
) -> tuple[dict, dict]:
    out_path = output_dir / f"{input_path.stem}.xyz"
    optimized_out_path = optimized_output_dir / f"{input_path.stem}.xyz"
    try:
        atoms = read(input_path, format="extxyz")
    except Exception as exc:
        return (
            _failed_row(input_path, out_path, None, exc),
            _optimized_summary_row(_failed_row(input_path, optimized_out_path, None, exc), None, None, None, None),
        )

    if max_atoms is not None and len(atoms) > max_atoms:
        error = f"crystal cell has {len(atoms)} atoms; cutoff is {max_atoms}"
        row = _force_summary_row(input_path, out_path, "skipped", len(atoms), None, error)
        optimized_row = _optimized_summary_row(
            _force_summary_row(input_path, optimized_out_path, "skipped", len(atoms), None, error),
            h_atoms=int(np.count_nonzero(atoms.symbols == "H")),
            optimization_converged=None,
            optimization_steps=None,
            max_h_force_ev_per_ang=None,
        )
        return row, optimized_row

    try:
        row = _calculate_and_write_forces(atoms.copy(), input_path, out_path)
    except Exception as exc:
        if out_path.exists():
            out_path.unlink()
        row = _failed_row(input_path, out_path, len(atoms), exc)
        optimized_row = _optimized_summary_row(
            _force_summary_row(
                input_path,
                optimized_out_path,
                "skipped",
                len(atoms),
                None,
                "unoptimized force calculation failed",
            ),
            h_atoms=int(np.count_nonzero(atoms.symbols == "H")),
            optimization_converged=None,
            optimization_steps=None,
            max_h_force_ev_per_ang=None,
        )
        return row, optimized_row

    try:
        optimized_row = _optimize_hydrogens_and_write_forces(
            atoms.copy(),
            input_path,
            optimized_out_path,
            h_fmax=h_fmax,
            h_steps=h_steps,
        )
    except Exception as exc:
        if optimized_out_path.exists():
            optimized_out_path.unlink()
        optimized_row = _optimized_summary_row(
            _failed_row(input_path, optimized_out_path, len(atoms), exc),
            h_atoms=int(np.count_nonzero(atoms.symbols == "H")),
            optimization_converged=False,
            optimization_steps=None,
            max_h_force_ev_per_ang=None,
        )

    return row, optimized_row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    parser.add_argument(
        "--optimized-output-dir",
        help="Directory for force-annotated structures after relaxing hydrogen positions only",
    )
    parser.add_argument("--limit", type=int, help="Only process the first N selected extxyz files")
    parser.add_argument("--prefix-filter", help="Only process extxyz files with ids starting with this prefix")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument(
        "--h-fmax",
        type=float,
        default=DEFAULT_H_FMAX,
        help="FIRE convergence threshold for hydrogen-only relaxation in eV/A",
    )
    parser.add_argument(
        "--h-steps",
        type=int,
        default=DEFAULT_H_STEPS,
        help="Maximum FIRE steps for hydrogen-only relaxation",
    )
    parser.add_argument(
        "--max-atoms",
        type=int,
        default=DEFAULT_MAX_ATOMS,
        help="Skip expanded unit cells above this atom count; use 0 to disable",
    )
    args = parser.parse_args()
    max_atoms = None if args.max_atoms == 0 else args.max_atoms

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    optimized_output_dir = (
        Path(args.optimized_output_dir)
        if args.optimized_output_dir
        else output_dir.with_name(f"{output_dir.name}_h_optimized")
    )
    if output_dir == optimized_output_dir:
        raise SystemExit("output_dir and optimized_output_dir must be different")

    output_dir.mkdir(parents=True, exist_ok=True)
    optimized_output_dir.mkdir(parents=True, exist_ok=True)
    _clear_previous_outputs(output_dir)
    _clear_previous_outputs(optimized_output_dir)

    input_files = sorted(input_dir.glob("*.xyz"))
    if args.prefix_filter:
        input_files = [path for path in input_files if path.stem.startswith(args.prefix_filter)]
    if args.limit is not None:
        input_files = input_files[: args.limit]
    if not input_files:
        raise SystemExit(f"No extxyz files found in {input_dir}")

    workers = max(1, min(args.workers, len(input_files)))
    rows = []
    optimized_rows = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                _compute_one,
                input_path,
                output_dir,
                optimized_output_dir,
                max_atoms,
                args.h_fmax,
                args.h_steps,
            ): input_path
            for input_path in input_files
        }
        for future in as_completed(futures):
            row, optimized_row = future.result()
            rows.append(row)
            optimized_rows.append(optimized_row)
            print(
                f"{row['status']}: {row['mol_id']} -> {row['extxyz_path']}; "
                f"h-optimized {optimized_row['status']} -> {optimized_row['extxyz_path']}"
            )

    summary = pd.DataFrame(rows).sort_values("mol_id")
    summary.to_csv(output_dir / SUMMARY_PATH, index=False)
    optimized_summary = pd.DataFrame(optimized_rows).sort_values("mol_id")
    optimized_summary.to_csv(optimized_output_dir / SUMMARY_PATH, index=False)
    if (summary["status"] == "failed").any() or (optimized_summary["status"] == "failed").any():
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
