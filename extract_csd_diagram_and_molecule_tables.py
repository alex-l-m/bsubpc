"""
extract_csd_diagram_and_molecule_tables.py

Goal
----
Dump TWO graph-like representations from the CSD Python API into CSV tables:

(1) The *internal 2D chemical diagram* graph stored on the entry:
    entry._entry.chemical_diagram_views().diagram()

    This is an undocumented SWIG-wrapped object (ChemistryLib.ChemicalDiagram).
    We learned by probing it that it provides:
      - natoms(), nbonds(), ncomponents()
      - atom(i).index(), atom(i).element(), atom(i).formal_charge(), atom(i).cyclic()
      - atom(i).centre_label()/side_label()/top_left_label()/top_right_label()
      - bond(i).atom1()/atom2() and bond(i).type().text()/value()

(2) The *crystal molecule* chemistry graph:
    entry.crystal.molecule

    This is the normal documented CCDC Molecule API.
    Important detail we learned from NAKMUI:
      - Some atoms can be "siteless" (atom.coordinates is None).
      - Those atoms are part of the chemical graph, but cannot appear in MOL2.

Outputs
-------
Writes 6 CSVs into out_dir (default: csd_tables/):

  csd_diagrams.csv
  csd_diagram_atoms.csv
  csd_diagram_bonds.csv

  csd_molecules.csv
  csd_molecule_atoms.csv   <-- now includes: has_coordinates, x, y, z, occupancy
  csd_molecule_bonds.csv

Usage
-----
  run_csd_python_api extract_csd_diagram_and_molecule_tables.py
  run_csd_python_api extract_csd_diagram_and_molecule_tables.py csd_info.csv csd_tables
"""

import os
import sys

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
    diagram_atom_rows = []
    diagram_bond_rows = []

    mol_rows = []
    mol_atom_rows = []
    mol_bond_rows = []

    for refcode in refcodes:
        entry = reader.entry(refcode)

        # -----------------------------
        # (A) Diagram tables
        # -----------------------------
        cdv = entry._entry.chemical_diagram_views()
        d = cdv.diagram()
        n_views = count_diagram_views(cdv)

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
        mol = entry.crystal.molecule

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

    # Write all 6 tables
    pd.DataFrame(diagram_rows).to_csv(os.path.join(outdir, "csd_diagrams.csv"), index=False)
    pd.DataFrame(diagram_atom_rows).to_csv(os.path.join(outdir, "csd_diagram_atoms.csv"), index=False)
    pd.DataFrame(diagram_bond_rows).to_csv(os.path.join(outdir, "csd_diagram_bonds.csv"), index=False)

    pd.DataFrame(mol_rows).to_csv(os.path.join(outdir, "csd_molecules.csv"), index=False)
    pd.DataFrame(mol_atom_rows).to_csv(os.path.join(outdir, "csd_molecule_atoms.csv"), index=False)
    pd.DataFrame(mol_bond_rows).to_csv(os.path.join(outdir, "csd_molecule_bonds.csv"), index=False)


if __name__ == "__main__":
    main()
