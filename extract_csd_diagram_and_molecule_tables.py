"""
extract_csd_diagram_and_molecule_tables.py

Goal
----
Dump TWO graph-like representations from the CSD Python API into CSV tables,
plus a per-component formula table that downstream filters use to gate mol
and cif outputs.

(1) The *internal 2D chemical diagram* graph stored on the entry:
    entry._entry.chemical_diagram_views().diagram()

    This is an undocumented SWIG-wrapped object (ChemistryLib.ChemicalDiagram).
    Probing it (see CLAUDE.md) revealed:
      - natoms(), nbonds(), ncomponents()
      - cdv.diagram_view(i) iterates *unique* component views (a unit cell
        with two identical molecules in the asym unit yields ONE view, not
        two — the diagram is deduplicated)
      - atom(i).index(), atom(i).element(), atom(i).formal_charge(), atom(i).cyclic()
      - bond(i).atom1()/atom2() and bond(i).type().text()/value()
      - n_formula_groups(), n_symbol_groups(), n_repeat_groups(),
        n_ligand_groups(); only `n_symbol_groups` ever fires on BsubPc data
        and it does *not* collapse atom counts.

(2) The *crystal molecule* chemistry graph:
    entry.crystal.molecule

    The normal documented CCDC Molecule API. `mol.components` decomposes it
    into connected components — each one a molecule in the asym unit.
    Important detail: some atoms can be "siteless" (atom.coordinates is
    None). Those atoms are part of the chemical graph but cannot appear in
    MOL2.

Outputs
-------
Writes 9 CSVs into out_dir (default: csd_tables/):

  csd_diagrams.csv
  csd_diagram_checks.csv         <-- formula-based, 'ok' iff every active
                                     molecule component has a matching
                                     diagram-view formula
  csd_diagram_atoms.csv
  csd_diagram_bonds.csv
  csd_diagram_view_formulas.csv  <-- new: csd_refcode, view_id, formula

  csd_molecules.csv
  csd_molecule_atoms.csv         <-- includes: has_coordinates, x, y, z, occupancy
  csd_molecule_bonds.csv
  csd_molecule_component_formulas.csv  <-- new: csd_refcode, component_id, formula

Usage
-----
  run_csd_python_api extract_csd_diagram_and_molecule_tables.py
  run_csd_python_api extract_csd_diagram_and_molecule_tables.py csd_info.csv csd_tables
"""

import os
import sys
from collections import Counter

import pandas as pd
import ccdc.io


def count_diagram_views(cdv):
    """Count cdv.diagram_view(i) views (i=0,1,2,...) until it errors."""
    i = 0
    while True:
        try:
            cdv.diagram_view(i)
        except Exception:
            return i
        i += 1


def _formula_string(elements) -> str:
    """Canonical sorted formula string, e.g. 'B1C24Cl13N6'."""
    counter = Counter(elements)
    return "".join(f"{el}{counter[el]}" for el in sorted(counter))


def _view_elements(view):
    """Iterate the element symbols of a diagram view's atoms."""
    for i in range(view.natoms()):
        a = view.atom(i)
        el = a.element()
        # `element()` is a ChemistryLib.Element; its `symbol()` accessor exists
        # in newer CCDC builds, but `str(el)` round-trips reliably as
        # "Symbol (Z)" so we strip the Z.
        if hasattr(el, "symbol") and callable(el.symbol):
            yield el.symbol()
        else:
            yield str(el).split(None, 1)[0]


def _component_formulas(components):
    """Sorted formula string per connected component (preserves order)."""
    return [_formula_string(a.atomic_symbol for a in c.atoms) for c in components]


def _view_formulas(cdv):
    """Sorted formula string per diagram view (== unique connected component)."""
    n_views = cdv.ncomponents()
    return [_formula_string(_view_elements(cdv.diagram_view(i))) for i in range(n_views)]


def _formula_check_status(component_formulas, view_formulas, has_formula_group):
    """Decide whether the active molecule is consistent with the diagram.

    Pass criterion: every active-molecule connected-component formula
    appears in at least one diagram view's formula. Diagram views are
    deduplicated, so multiple identical components all pass against the same
    view. We do not check the converse direction (extra views in the
    diagram), since the chemist may have drawn a disorder alternate or
    co-crystal solvent that is absent from the selected active state.

    If any view uses a formula group (collapsed atom shorthand like '-CF3'),
    the diagram-side formula is incomplete and the check is unreliable; we
    return ``'unknown'`` and let downstream filters treat that as a pass.
    """
    if has_formula_group:
        return "unknown", "diagram contains a formula group; atom count incomplete"
    missing = [
        formula for formula in component_formulas if formula not in view_formulas
    ]
    if missing:
        # Distinct unmatched formulas, comma-joined.
        unmatched = sorted(set(missing))
        return (
            "failed",
            f"active component formula(s) missing from diagram: {', '.join(unmatched)}",
        )
    return "ok", ""


def main():
    inpath = sys.argv[1] if len(sys.argv) > 1 else "csd_info.csv"
    outdir = sys.argv[2] if len(sys.argv) > 2 else "csd_tables"
    os.makedirs(outdir, exist_ok=True)

    intbl = pd.read_csv(inpath)

    # Remove rows with missing "structure_entry" column
    intbl = intbl[intbl["structure_entry"].notna()]

    # Pull refcodes from column "structure_entry"
    refcodes = intbl["structure_entry"].tolist()

    reader = ccdc.io.EntryReader("CSD")

    diagram_rows = []
    diagram_check_rows = []
    diagram_atom_rows = []
    diagram_bond_rows = []
    diagram_view_formula_rows = []

    mol_rows = []
    mol_atom_rows = []
    mol_bond_rows = []
    mol_component_formula_rows = []

    for refcode in refcodes:
        entry = reader.entry(refcode)

        # -----------------------------
        # (A) Diagram tables
        # -----------------------------
        cdv = entry._entry.chemical_diagram_views()
        d = cdv.diagram()
        n_views = count_diagram_views(cdv)
        active_mol = entry.crystal.molecule
        disordered_mol = entry.crystal.disordered_molecule

        view_formulas = _view_formulas(cdv)
        active_components = list(active_mol.components)
        component_formulas = _component_formulas(active_components)
        for view_id, formula in enumerate(view_formulas):
            view = cdv.diagram_view(view_id)
            diagram_view_formula_rows.append(
                {
                    "csd_refcode": refcode,
                    "view_id": view_id,
                    "formula": formula,
                    "natoms": view.natoms(),
                    "n_formula_groups": view.n_formula_groups(),
                    "n_symbol_groups": view.n_symbol_groups(),
                    "n_ligand_groups": view.n_ligand_groups(),
                    "n_repeat_groups": view.n_repeat_groups(),
                }
            )
        for component_id, formula in enumerate(component_formulas):
            mol_component_formula_rows.append(
                {
                    "csd_refcode": refcode,
                    "component_id": component_id,
                    "formula": formula,
                    "natoms": sum(1 for _ in active_components[component_id].atoms),
                }
            )

        has_formula_group = any(
            cdv.diagram_view(i).n_formula_groups() > 0 for i in range(n_views)
        )
        check_status, check_reason = _formula_check_status(
            component_formulas, view_formulas, has_formula_group
        )

        diagram_rows.append(
            {
                "csd_refcode": refcode,
                "diagram_class": type(d).__name__,
                "diagram_natoms": d.natoms(),
                "diagram_nbonds": d.nbonds(),
                "diagram_ncomponents": d.ncomponents(),
                "diagram_n_ligand_groups": d.n_ligand_groups(),
                "diagram_n_symbol_groups": d.n_symbol_groups(),
                "diagram_n_repeat_groups": d.n_repeat_groups(),
                "diagram_n_formula_groups": d.n_formula_groups(),
                "diagram_n_views": n_views,
            }
        )
        # The "ok" gate now is per-component formula match: every connected
        # component of the active molecule must have a diagram view with the
        # same element-count signature. This still uses crystal.molecule (the
        # active disorder state) — search.py writes that selection back into
        # the crystal before emitting the cif, so it is the structure
        # downstream pipelines actually consume.
        diagram_check_rows.append(
            {
                "csd_refcode": refcode,
                "status": check_status,
                "reason": check_reason,
                "diagram_natoms": d.natoms(),
                "diagram_nbonds": d.nbonds(),
                "active_molecule_natoms": len(active_mol.atoms),
                "active_molecule_nbonds": len(active_mol.bonds),
                "disordered_molecule_natoms": len(disordered_mol.atoms),
                "disordered_molecule_nbonds": len(disordered_mol.bonds),
                "active_component_formulas": "|".join(component_formulas),
                "diagram_view_formulas": "|".join(view_formulas),
            }
        )

        for i in range(d.natoms()):
            a = d.atom(i)
            diagram_atom_rows.append(
                {
                    "csd_refcode": refcode,
                    "diagram_atom_id": a.index(),
                    "element": a.element(),
                    "formal_charge": a.formal_charge(),
                    "cyclic": a.cyclic(),
                    "centre_label": str(a.centre_label()),
                    "side_label": str(a.side_label()),
                    "top_left_label": str(a.top_left_label()),
                    "top_right_label": str(a.top_right_label()),
                }
            )

        for i in range(d.nbonds()):
            b = d.bond(i)
            a1 = b.atom1()
            a2 = b.atom2()
            bt = b.type()
            diagram_bond_rows.append(
                {
                    "csd_refcode": refcode,
                    "diagram_bond_id": i,
                    "atom1_id": a1.index(),
                    "atom2_id": a2.index(),
                    "bond_type_text": bt.text(),
                    "bond_type_value": bt.value(),
                }
            )

        # -----------------------------
        # (B) Molecule tables
        # -----------------------------
        # Per your preference: crystal.molecule is fine; key is NOT to use disordered_molecule.
        mol = active_mol

        mol_rows.append(
            {
                "csd_refcode": refcode,
                "molecule_class": type(mol).__name__,
                "mol_natoms": len(mol.atoms),
                "mol_nbonds": len(mol.bonds),
                "mol_formula": getattr(mol, "formula", None),
                "mol_formal_charge": getattr(mol, "formal_charge", None),
                "mol_smiles": getattr(mol, "smiles", None),
            }
        )

        for a in mol.atoms:
            coords = a.coordinates  # can be None (siteless atom)
            has_coords = int(coords is not None)
            mol_atom_rows.append(
                {
                    "csd_refcode": refcode,
                    "mol_atom_id": a.index,
                    "element": a.atomic_symbol,
                    "formal_charge": a.formal_charge,
                    "label": a.label,
                    # New columns: everything needed to detect "siteless" atoms later.
                    "has_coordinates": has_coords,
                    "x": None if coords is None else coords.x,
                    "y": None if coords is None else coords.y,
                    "z": None if coords is None else coords.z,
                    "occupancy": getattr(a, "occupancy", None),
                }
            )

        for i, b in enumerate(mol.bonds):
            a1, a2 = b.atoms

            # For molecule bonds, the simplest reliable interface is:
            #   str(b.bond_type)  -> "Single", "Aromatic", ...
            # and for numeric code we can use the underlying _bond_type if available.
            bt = b.bond_type
            bt_text = str(bt)
            bt_val = bt._bond_type.value() if hasattr(bt, "_bond_type") else None

            mol_bond_rows.append(
                {
                    "csd_refcode": refcode,
                    "mol_bond_id": i,
                    "atom1_id": a1.index,
                    "atom2_id": a2.index,
                    "bond_type_text": bt_text,
                    "bond_type_value": bt_val,
                }
            )

    # Write all output tables
    pd.DataFrame(diagram_rows).to_csv(os.path.join(outdir, "csd_diagrams.csv"), index=False)
    pd.DataFrame(diagram_check_rows).to_csv(os.path.join(outdir, "csd_diagram_checks.csv"), index=False)
    pd.DataFrame(diagram_atom_rows).to_csv(os.path.join(outdir, "csd_diagram_atoms.csv"), index=False)
    pd.DataFrame(diagram_bond_rows).to_csv(os.path.join(outdir, "csd_diagram_bonds.csv"), index=False)
    pd.DataFrame(diagram_view_formula_rows).to_csv(
        os.path.join(outdir, "csd_diagram_view_formulas.csv"), index=False
    )

    pd.DataFrame(mol_rows).to_csv(os.path.join(outdir, "csd_molecules.csv"), index=False)
    pd.DataFrame(mol_atom_rows).to_csv(os.path.join(outdir, "csd_molecule_atoms.csv"), index=False)
    pd.DataFrame(mol_bond_rows).to_csv(os.path.join(outdir, "csd_molecule_bonds.csv"), index=False)
    pd.DataFrame(mol_component_formula_rows).to_csv(
        os.path.join(outdir, "csd_molecule_component_formulas.csv"), index=False
    )


if __name__ == "__main__":
    main()
