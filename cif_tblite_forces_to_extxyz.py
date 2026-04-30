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
from ase.io import read, write
from tblite.ase import TBLite


DEFAULT_WORKERS = os.cpu_count() or 1
SUMMARY_PATH = "tblite_force_summary.csv"
DEFAULT_MAX_ATOMS = 300


def _clear_previous_outputs(output_dir: Path) -> None:
    for path in output_dir.glob("*.xyz"):
        path.unlink()
    summary_path = output_dir / SUMMARY_PATH
    if summary_path.exists():
        summary_path.unlink()


def _compute_one(input_path: Path, output_dir: Path, max_atoms: int | None) -> dict:
    out_path = output_dir / f"{input_path.stem}.xyz"
    try:
        atoms = read(input_path, format="extxyz")
        if max_atoms is not None and len(atoms) > max_atoms:
            return {
                "mol_id": input_path.stem,
                "input_path": str(input_path),
                "extxyz_path": str(out_path),
                "status": "skipped",
                "n_atoms": len(atoms),
                "max_force_ev_per_ang": None,
                "error": f"crystal cell has {len(atoms)} atoms; cutoff is {max_atoms}",
            }

        atoms.calc = TBLite(
            method="GFN2-xTB",
            charge=0,
            multiplicity=1,
            verbosity=0,
        )

        forces = atoms.get_forces()

        write(
            out_path,
            atoms,
            format="extxyz",
            columns=["symbols", "positions", "forces"],
            write_results=True,
        )
        return {
            "mol_id": input_path.stem,
            "input_path": str(input_path),
            "extxyz_path": str(out_path),
            "status": "ok",
            "n_atoms": len(atoms),
            "max_force_ev_per_ang": float(np.linalg.norm(forces, axis=1).max()),
            "error": "",
        }
    except Exception as exc:
        if out_path.exists():
            out_path.unlink()
        return {
            "mol_id": input_path.stem,
            "input_path": str(input_path),
            "extxyz_path": str(out_path),
            "status": "failed",
            "n_atoms": None,
            "max_force_ev_per_ang": None,
            "error": f"{type(exc).__name__}: {exc}",
        }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    parser.add_argument("--limit", type=int, help="Only process the first N selected extxyz files")
    parser.add_argument("--prefix-filter", help="Only process extxyz files with ids starting with this prefix")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
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
    output_dir.mkdir(parents=True, exist_ok=True)
    _clear_previous_outputs(output_dir)

    input_files = sorted(input_dir.glob("*.xyz"))
    if args.prefix_filter:
        input_files = [path for path in input_files if path.stem.startswith(args.prefix_filter)]
    if args.limit is not None:
        input_files = input_files[: args.limit]
    if not input_files:
        raise SystemExit(f"No extxyz files found in {input_dir}")

    workers = max(1, min(args.workers, len(input_files)))
    rows = []
    with ProcessPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(_compute_one, input_path, output_dir, max_atoms): input_path
            for input_path in input_files
        }
        for future in as_completed(futures):
            row = future.result()
            rows.append(row)
            print(f"{row['status']}: {row['mol_id']} -> {row['extxyz_path']}")

    summary = pd.DataFrame(rows).sort_values("mol_id")
    summary.to_csv(output_dir / SUMMARY_PATH, index=False)
    if (summary["status"] == "failed").any():
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
