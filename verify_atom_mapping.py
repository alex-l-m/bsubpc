#!/usr/bin/env python3
"""Verify molecule-to-CIF atom mappings by comparing pairwise distances."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from ase.io.cif import parse_cif
from rdkit.Chem import MolFromMol2File, MolFromMolFile


MAPPING_PATH = Path("atom_mapping.csv")
MOL_DIR = Path("crystal_geom")
FALLBACK_MOL2_DIR = Path("csd_molecules")
CIF_DIR = Path("csd_crystals")
OUTPUT_PATH = Path("atom_mapping_check.csv")
TOLERANCE = 1e-3


def _rdkit_positions(mol_path: Path) -> tuple[list[str], np.ndarray]:
    if mol_path.suffix.lower() == ".mol2":
        mol = MolFromMol2File(str(mol_path), removeHs=False, sanitize=False)
    else:
        mol = MolFromMolFile(str(mol_path), removeHs=False, sanitize=False)
    if mol is None:
        raise ValueError(f"RDKit could not read {mol_path}")
    conformer = mol.GetConformer()
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    positions = np.array(
        [
            [
                conformer.GetAtomPosition(atom_index).x,
                conformer.GetAtomPosition(atom_index).y,
                conformer.GetAtomPosition(atom_index).z,
            ]
            for atom_index in range(mol.GetNumAtoms())
        ],
        dtype=float,
    )
    return symbols, positions


def _cif_site_positions(cif_path: Path) -> tuple[list[str], np.ndarray]:
    blocks = [block for block in parse_cif(str(cif_path)) if block.has_structure()]
    if len(blocks) != 1:
        raise ValueError(f"expected exactly one structural CIF block in {cif_path}, found {len(blocks)}")
    atoms = blocks[0].get_unsymmetrized_structure()
    return atoms.get_chemical_symbols(), atoms.get_positions()


def _distance_matrix(positions: np.ndarray) -> np.ndarray:
    delta = positions[:, None, :] - positions[None, :, :]
    return np.linalg.norm(delta, axis=2)


def _molecule_path(mol_id: str) -> Path:
    mol_path = MOL_DIR / f"{mol_id}.mol"
    if mol_path.exists():
        return mol_path
    mol2_path = FALLBACK_MOL2_DIR / f"{mol_id}.mol2"
    if mol2_path.exists():
        return mol2_path
    raise FileNotFoundError(f"no molecule file found for {mol_id}")


def check_one(mol_id: str, rows: pd.DataFrame) -> dict:
    mol_path = _molecule_path(mol_id)
    cif_path = CIF_DIR / f"{mol_id}.cif"
    if not cif_path.exists():
        raise FileNotFoundError(f"no CIF file found for {mol_id}: {cif_path}")

    mol_symbols, mol_positions = _rdkit_positions(mol_path)
    cif_symbols, cif_positions = _cif_site_positions(cif_path)

    rows = rows.sort_values(["molecule_atom_idx", "crystal_atom_idx"])
    mol_idx = rows["molecule_atom_idx"].to_numpy(dtype=int)
    cif_idx = rows["crystal_atom_idx"].to_numpy(dtype=int)

    if np.any(mol_idx < 0) or np.any(mol_idx >= len(mol_symbols)):
        raise ValueError("mapping contains molecule_atom_idx outside the molecule atom range")
    if np.any(cif_idx < 0) or np.any(cif_idx >= len(cif_symbols)):
        raise ValueError("mapping contains crystal_atom_idx outside the CIF atom-site range")

    element_mismatches = [
        f"{int(m)}:{mol_symbols[int(m)]}!={cif_symbols[int(c)]}"
        for m, c in zip(mol_idx, cif_idx, strict=True)
        if mol_symbols[int(m)] != cif_symbols[int(c)]
    ]
    if element_mismatches:
        raise ValueError(f"mapped elements differ: {', '.join(element_mismatches[:10])}")

    mol_distances = _distance_matrix(mol_positions[mol_idx])
    cif_distances = _distance_matrix(cif_positions[cif_idx])
    absolute_difference = np.abs(mol_distances - cif_distances)
    max_abs_distance_error = float(np.max(absolute_difference)) if len(rows) else 0.0

    return {
        "mol_id": mol_id,
        "status": "ok" if max_abs_distance_error <= TOLERANCE else "failed",
        "n_mapped_atoms": len(rows),
        "molecule_path": str(mol_path),
        "cif_path": str(cif_path),
        "max_abs_distance_error": max_abs_distance_error,
        "tolerance": TOLERANCE,
        "error": "" if max_abs_distance_error <= TOLERANCE else "distance matrices differ",
    }


def main() -> int:
    mapping = pd.read_csv(MAPPING_PATH)
    required_columns = {"mol_id", "molecule_atom_idx", "crystal_atom_idx"}
    missing = required_columns.difference(mapping.columns)
    if missing:
        raise ValueError(f"{MAPPING_PATH} is missing required columns: {', '.join(sorted(missing))}")

    report_rows = []
    for mol_id, rows in mapping.groupby("mol_id", sort=True):
        try:
            report_rows.append(check_one(str(mol_id), rows))
        except Exception as exc:
            report_rows.append(
                {
                    "mol_id": mol_id,
                    "status": "failed",
                    "n_mapped_atoms": len(rows),
                    "molecule_path": "",
                    "cif_path": "",
                    "max_abs_distance_error": np.nan,
                    "tolerance": TOLERANCE,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )

    report = pd.DataFrame(report_rows)
    report.to_csv(OUTPUT_PATH, index=False)
    if not report.empty and (report["status"] != "ok").any():
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
