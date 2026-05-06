#!/usr/bin/env python3
"""Process BsubPc search results into labeled, canonicalized outputs.

This script:
- reads ``search_results.csv``;
- when a ``.mol2`` is available and the matched component's element-count
  signature matches some CSD chemical-diagram view's formula, sanitizes /
  canonicalizes it and writes a ``crystal_geom/<id>.mol`` for downstream
  sTDA work;
- when a ``.cif`` is available and the active molecule's per-component
  formulas all appear among the CSD diagram views, loads the P1-expanded
  crystal cell and writes a ``crystal_cells/<id>.xyz`` for downstream force
  work;
- writes one diagram svg per successfully canonicalized molecule.

Mol-output and crystal-cell-output succeed/fail independently: a CSD entry
with a problematic unit cell can still contribute its mol2 to sTDA, and an
entry whose mol2 cannot be sanitized can still contribute its cif to force
analysis. Both gates are formula-based — the matched mol2 component must
have a diagram-view twin, and every active-mol component must have one for
the cif. If a diagram view uses a formula group (collapsed atom shorthand)
the diagram-side atom count is incomplete; in that case the gate is
disabled for that entry and the mol/cif are written without the consistency
check. Files that fail any per-step processing are skipped and logged in
``processing_summary.csv``.
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
DIAGRAM_VIEW_FORMULAS_PATH = Path("csd_tables/csd_diagram_view_formulas.csv")
STRUCTURE_DIR = Path("crystal_geom")
CRYSTAL_DIR = Path("crystal_cells")
DIAGRAM_DIR = Path("diagrams")
ATOM_TABLE_PATH = Path("atom_table.csv")
SUMMARY_PATH = Path("processing_summary.csv")


def clear_previous_outputs() -> None:
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


def load_diagram_view_formulas(path: Path) -> dict[str, set[str]]:
    """Return ``refcode -> {view_formula, …}`` from the diagram extractor output.

    Entries where any view used a formula group (collapsed atom shorthand)
    cannot be reliably formula-matched against the mol2; they are mapped to
    ``None`` so the caller can skip the check rather than reject everything.
    """
    if not path.exists():
        return {}
    table = pd.read_csv(path)
    grouped: dict[str, set[str] | None] = {}
    for refcode, sub in table.groupby("csd_refcode"):
        if (sub["n_formula_groups"] > 0).any():
            grouped[refcode] = None
        else:
            grouped[refcode] = set(sub["formula"].astype(str))
    return grouped


def passing_csd_diagram_ids(path: Path) -> set[str]:
    """Refcodes whose active molecule passes the per-component formula check.

    ``status == 'unknown'`` means a diagram formula group prevented the check;
    we treat that as a pass to avoid losing entries we cannot validate.
    """
    if not path.exists():
        return set()
    checks = pd.read_csv(path)
    return set(checks.loc[checks["status"].isin(["ok", "unknown"]), "csd_refcode"])


def load_crystal(path: Path):
    """Load the selected P1 periodic unit cell used for force calculation."""
    atoms = read(path, format="cif")
    atoms.set_pbc((True, True, True))
    return atoms


def process_mol2(mol2_path: Path) -> tuple[Mol, str]:
    """Canonicalize a mol2 and write the per-molecule mol + diagram svg."""
    mol = load_molecule(mol2_path)
    mol_id = mol.GetProp("_Name")
    MolToMolFile(mol, str(STRUCTURE_DIR / f"{mol_id}.mol"))
    diagram_error = ""
    try:
        write_diagram(mol, str(DIAGRAM_DIR / f"{mol_id}.svg"))
    except Exception as exc:
        diagram_error = f"{type(exc).__name__}: {exc}"
    return mol, diagram_error


def process_cif(cif_path: Path, mol_id: str) -> None:
    """Read the P1 cif and write a ``crystal_cells/<id>.xyz`` extxyz file."""
    atoms = load_crystal(cif_path)
    write(CRYSTAL_DIR / f"{mol_id}.xyz", atoms, format="extxyz")


def remove_mol_outputs(mol_id: str) -> None:
    for path in (
        STRUCTURE_DIR / f"{mol_id}.mol",
        DIAGRAM_DIR / f"{mol_id}.svg",
    ):
        path.unlink(missing_ok=True)


def remove_crystal_outputs(mol_id: str) -> None:
    (CRYSTAL_DIR / f"{mol_id}.xyz").unlink(missing_ok=True)


def _format_error(exc: BaseException) -> str:
    return f"{type(exc).__name__}: {exc}"


def main():
    clear_previous_outputs()
    search_results = pd.read_csv(CSV_PATH)
    passing_diagram_ids = passing_csd_diagram_ids(DIAGRAM_CHECKS_PATH)
    diagram_view_formulas = load_diagram_view_formulas(DIAGRAM_VIEW_FORMULAS_PATH)

    report_rows = []
    mols = []
    for row in search_results.itertuples(index=False):
        mol_id = row.mol_id
        mol2_path_str = row.mol2_path if isinstance(row.mol2_path, str) else ""
        cif_path_str = row.cif_path if isinstance(row.cif_path, str) else ""
        mol2_formula = row.mol2_formula if isinstance(row.mol2_formula, str) else ""

        mol_status = "skipped"
        mol_error = ""
        diagram_error = ""
        if mol2_path_str:
            views = diagram_view_formulas.get(mol_id)
            if views is not None and mol2_formula and mol2_formula not in views:
                mol_status = "failed"
                mol_error = (
                    f"matched-component formula {mol2_formula!r} not present "
                    f"in CSD diagram views {sorted(views)!r}"
                )
                remove_mol_outputs(mol_id)
            else:
                mol2_path = Path(mol2_path_str)
                try:
                    mol, diagram_error = process_mol2(mol2_path)
                    mol_status = "ok"
                    mols.append(mol)
                except Exception as exc:
                    mol_status = "failed"
                    mol_error = _format_error(exc)
                    remove_mol_outputs(mol_id)

        crystal_status = "skipped"
        crystal_error = ""
        if cif_path_str:
            if mol_id not in passing_diagram_ids:
                crystal_status = "failed"
                crystal_error = (
                    "active-molecule component formulas missing from CSD diagram views"
                )
                remove_crystal_outputs(mol_id)
            else:
                cif_path = Path(cif_path_str)
                try:
                    process_cif(cif_path, mol_id)
                    crystal_status = "ok"
                except Exception as exc:
                    crystal_status = "failed"
                    crystal_error = _format_error(exc)
                    remove_crystal_outputs(mol_id)

        report_rows.append(
            {
                "mol_id": mol_id,
                "mol2_path": mol2_path_str,
                "cif_path": cif_path_str,
                "mol2_formula": mol2_formula,
                "mol_status": mol_status,
                "mol_error": mol_error,
                "diagram_error": diagram_error,
                "crystal_status": crystal_status,
                "crystal_error": crystal_error,
                "crystal_path": (
                    str(CRYSTAL_DIR / f"{mol_id}.xyz") if crystal_status == "ok" else ""
                ),
            }
        )
    if mols:
        mols_frame = mols2tbl(mols)
    else:
        mols_frame = pd.DataFrame()
    mols_frame.to_csv(ATOM_TABLE_PATH, index=False)
    report_frame = pd.DataFrame(report_rows)
    report_frame.to_csv(SUMMARY_PATH, index=False)


if __name__ == "__main__":
    raise SystemExit(main())
