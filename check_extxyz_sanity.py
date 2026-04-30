#!/usr/bin/env python3
"""Check force-annotated extxyz files for unreasonable positions or forces."""

from pathlib import Path

import numpy as np
import pandas as pd
from ase.io import read


INPUT_DIR = Path("crystal_forces")
OUTPUT_PATH = Path("extxyz_sanity_check.csv")
MAX_ABS_POSITION_ANG = 1.0e4
MAX_FORCE_EV_PER_ANG = 100.0
MIN_DISTANCE_ANG = 0.5


def check_one(path: Path) -> dict:
    try:
        atoms = read(path, format="extxyz")
        positions = atoms.get_positions()
        forces = atoms.get_forces()
        force_norms = np.linalg.norm(forces, axis=1)

        finite_positions = bool(np.isfinite(positions).all())
        finite_forces = bool(np.isfinite(forces).all())
        max_abs_position = float(np.max(np.abs(positions))) if len(atoms) else 0.0
        max_force = float(np.max(force_norms)) if len(atoms) else 0.0
        min_distance = np.nan
        if len(atoms) > 1:
            distances = atoms.get_all_distances(mic=True)
            np.fill_diagonal(distances, np.nan)
            min_distance = float(np.nanmin(distances))

        problems = []
        if len(atoms) == 0:
            problems.append("no atoms")
        if not finite_positions:
            problems.append("nonfinite positions")
        if not finite_forces:
            problems.append("nonfinite forces")
        if max_abs_position > MAX_ABS_POSITION_ANG:
            problems.append(f"max abs position exceeds {MAX_ABS_POSITION_ANG:g} A")
        if max_force > MAX_FORCE_EV_PER_ANG:
            problems.append(f"max force exceeds {MAX_FORCE_EV_PER_ANG:g} eV/A")
        if np.isfinite(min_distance) and min_distance < MIN_DISTANCE_ANG:
            problems.append(f"min distance below {MIN_DISTANCE_ANG:g} A")

        return {
            "mol_id": path.stem,
            "path": str(path),
            "status": "ok" if not problems else "failed",
            "n_atoms": len(atoms),
            "max_abs_position_ang": max_abs_position,
            "max_force_ev_per_ang": max_force,
            "min_distance_ang": min_distance,
            "error": "; ".join(problems),
        }
    except Exception as exc:
        return {
            "mol_id": path.stem,
            "path": str(path),
            "status": "failed",
            "n_atoms": None,
            "max_abs_position_ang": np.nan,
            "max_force_ev_per_ang": np.nan,
            "min_distance_ang": np.nan,
            "error": f"{type(exc).__name__}: {exc}",
        }


def main() -> int:
    paths = sorted(INPUT_DIR.glob("*.xyz"))
    if not paths:
        raise SystemExit(f"No extxyz files found in {INPUT_DIR}")

    report = pd.DataFrame(check_one(path) for path in paths)
    report.to_csv(OUTPUT_PATH, index=False)
    if (report["status"] != "ok").any():
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
