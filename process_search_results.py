#!/usr/bin/env python3
"""Process BsubPc search results into labeled, canonicalized outputs.

This script:
- reads ``search_results.csv``;
- loads each unique ``.mol2`` listed there;
- writes an SVG diagram
- appends one row per atom to ``atom_table.csv`
- saves a canonicalized ``.mol`` file for each successfully processed structure.

Files that RDKit cannot read/sanitize, or that fail any later processing step, are
skipped and logged in ``processing_summary.csv``.
"""

from pathlib import Path

import pandas as pd
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToMolFile, MolFromMol2File

from atom_table import mols2tbl
from canonicalize import canonicalize_bsubpc
from diagrams import write_diagram


CSV_PATH = Path("search_results.csv")
STRUCTURE_DIR = Path("crystal_geom")
DIAGRAM_DIR = Path("diagrams")
ATOM_TABLE_PATH = Path("atom_table.csv")
SUMMARY_PATH = Path("processing_summary.csv")


def unique_mol2_paths(csv_path: Path) -> list[Path]:
    """Return unique mol2 paths from the CSV, preserving first-seen order."""
    frame = pd.read_csv(csv_path)
    if "mol2_path" not in frame.columns:
        raise ValueError(f"{csv_path} is missing required column 'mol2_path'")

    seen: set[Path] = set()
    unique_paths: list[Path] = []
    for raw_path in frame["mol2_path"]:
        if pd.isna(raw_path):
            continue
        rel_path = Path(str(raw_path))
        if rel_path not in seen:
            seen.add(rel_path)
            unique_paths.append(rel_path)
    return unique_paths


def clear_previous_outputs() -> None:
    """Remove stale outputs so each run reflects only the current processing."""
    STRUCTURE_DIR.mkdir(exist_ok=True)
    DIAGRAM_DIR.mkdir(exist_ok=True)

    for path in STRUCTURE_DIR.glob("*.mol"):
        path.unlink()
    for path in DIAGRAM_DIR.glob("*.svg"):
        path.unlink()
    for path in (ATOM_TABLE_PATH, SUMMARY_PATH):
        if path.exists():
            path.unlink()


def load_molecule(path: Path):
    """Load a mol2 file and require RDKit sanitization to succeed."""
    mol = MolFromMol2File(str(path), removeHs=False, sanitize=True)
    if mol is None:
        raise ValueError("RDKit could not read/sanitize file")
    mol.SetProp("_Name", path.stem)
    canonicalize_bsubpc(mol)
    return mol


def process_one(in_path: Path) -> Mol:
    """Process a single structure"""
    mol = load_molecule(in_path)

    mol_id = mol.GetProp("_Name")
    MolToMolFile(mol, str(STRUCTURE_DIR / f"{mol_id}.mol"))
    write_diagram(mol, str(DIAGRAM_DIR / f"{mol_id}.svg"))
    return mol


def main():
    clear_previous_outputs()
    rel_paths = unique_mol2_paths(CSV_PATH)

    report_rows = []
    mols = []
    for rel_path in rel_paths:
        mol_id = Path(rel_path).stem
        try:
            mol = process_one(rel_path)
            mol_id = mol.GetProp("_Name")
            report_rows.append(
                {
                    "mol_id": mol_id,
                    "mol2_path": str(rel_path),
                    "status": "ok",
                    "error": "",
                }
            )
            mols.append(mol)
        except Exception as exc:
            report_rows.append(
                {
                    "mol_id": mol_id,
                    "mol2_path": str(rel_path),
                    "status": "skipped",
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
    mols_frame = mols2tbl(mols)
    mols_frame.to_csv(ATOM_TABLE_PATH, index=False)
    report_frame = pd.DataFrame(report_rows)
    report_frame.to_csv(SUMMARY_PATH, index=False)


if __name__ == "__main__":
    raise SystemExit(main())
