#!/usr/bin/env python3
"""Process BsubPc search results into labeled, canonicalized outputs.

This script:
- reads ``search_results.csv``;
- loads each ``.mol2`` listed there;
- filters structures whose CSD diagram atom count does not match the full structure;
- writes an SVG diagram;
- appends one row per atom to ``atom_table.csv`
- saves a canonicalized ``.mol`` file for each successfully processed structure.
- saves a symmetry-expanded crystal extxyz file for each successfully processed structure.

Files that RDKit cannot read/sanitize, or that fail any later processing step, are
skipped and logged in ``processing_summary.csv``.
"""

from pathlib import Path

import pandas as pd
from ase.io import read, write
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolToMolFile, MolFromMol2File

from atom_table import mols2tbl
from canonicalize import canonicalize_bsubpc
from diagrams import write_diagram


CSV_PATH = Path("search_results.csv")
DIAGRAM_CHECKS_PATH = Path("csd_tables/csd_diagram_checks.csv")
STRUCTURE_DIR = Path("crystal_geom")
CRYSTAL_DIR = Path("crystal_cells")
DIAGRAM_DIR = Path("diagrams")
ATOM_TABLE_PATH = Path("atom_table.csv")
SUMMARY_PATH = Path("processing_summary.csv")


def clear_previous_outputs() -> None:
    """Remove stale outputs so each run reflects only the current processing."""
    STRUCTURE_DIR.mkdir(exist_ok=True)
    CRYSTAL_DIR.mkdir(exist_ok=True)
    DIAGRAM_DIR.mkdir(exist_ok=True)

    for path in STRUCTURE_DIR.glob("*.mol"):
        path.unlink()
    for path in CRYSTAL_DIR.glob("*.xyz"):
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


def passing_csd_diagram_ids(path: Path) -> set[str]:
    """Return IDs whose CSD diagram atom count matches the full structure."""
    checks = pd.read_csv(path)
    return set(checks.loc[checks["status"] == "ok", "csd_refcode"])


def load_crystal(path: Path):
    """Load the symmetry-expanded periodic unit cell used for force calculation."""
    atoms = read(path, format="cif")
    atoms.set_pbc((True, True, True))
    return atoms


def process_one(mol2_path: Path, cif_path: Path) -> tuple[Mol, str]:
    """Process one molecule and its matching crystal cell."""
    mol = load_molecule(mol2_path)
    atoms = load_crystal(cif_path)

    mol_id = mol.GetProp("_Name")
    MolToMolFile(mol, str(STRUCTURE_DIR / f"{mol_id}.mol"))
    write(CRYSTAL_DIR / f"{mol_id}.xyz", atoms, format="extxyz")
    diagram_error = ""
    try:
        write_diagram(mol, str(DIAGRAM_DIR / f"{mol_id}.svg"))
    except Exception as exc:
        diagram_error = f"{type(exc).__name__}: {exc}"
    return mol, diagram_error


def remove_structure_outputs(mol_id: str) -> None:
    """Remove partial paired outputs for a skipped structure."""
    for path in (
        STRUCTURE_DIR / f"{mol_id}.mol",
        CRYSTAL_DIR / f"{mol_id}.xyz",
        DIAGRAM_DIR / f"{mol_id}.svg",
    ):
        path.unlink(missing_ok=True)


def main():
    clear_previous_outputs()
    # search.py writes search_results.csv with a non-null, unique mol2_path
    # for each accepted CSD refcode.
    search_results = pd.read_csv(CSV_PATH)
    passing_diagram_ids = passing_csd_diagram_ids(DIAGRAM_CHECKS_PATH)

    report_rows = []
    mols = []
    for row in search_results.itertuples(index=False):
        rel_path = Path(row.mol2_path)
        cif_path = Path(row.cif_path)
        mol_id = row.mol_id
        try:
            if mol_id not in passing_diagram_ids:
                raise ValueError("CSD diagram atom count does not match full structure atom count")
            mol, diagram_error = process_one(rel_path, cif_path)
            mol_id = mol.GetProp("_Name")
            report_rows.append(
                {
                    "mol_id": mol_id,
                    "mol2_path": str(rel_path),
                    "crystal_path": str(CRYSTAL_DIR / f"{mol_id}.xyz"),
                    "status": "ok",
                    "error": "",
                    "diagram_error": diagram_error,
                }
            )
            mols.append(mol)
        except Exception as exc:
            remove_structure_outputs(mol_id)
            report_rows.append(
                {
                    "mol_id": mol_id,
                    "mol2_path": str(rel_path),
                    "crystal_path": "",
                    "status": "skipped",
                    "error": f"{type(exc).__name__}: {exc}",
                    "diagram_error": "",
                }
            )
    mols_frame = mols2tbl(mols)
    mols_frame.to_csv(ATOM_TABLE_PATH, index=False)
    report_frame = pd.DataFrame(report_rows)
    report_frame.to_csv(SUMMARY_PATH, index=False)


if __name__ == "__main__":
    raise SystemExit(main())
