from pathlib import Path
from rdkit.Chem import rdDepictor
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolops import FragmentOnBonds, RemoveHs, AddHs

from match import assert_labeled


def split_subpc_at_apical_bond(mol: Mol) -> Mol:
    """Fragment a labeled BsubPc at its boron-axial bond and return the
    disconnected (axial-substituent, bowl) Mol with attachment dummies
    at the cleavage."""
    assert_labeled(mol)
    boron_idx = next(
        a.GetIdx() for a in mol.GetAtoms()
        if a.HasProp("bsubpc_label") and a.GetProp("bsubpc_label") == "boron"
    )
    # The boron's three template-labeled neighbors are the pyrrole
    # nitrogens; the only neighbor without `bsubpc_label` is the axial
    # substituent atom.
    axial_idx = next(
        nbr.GetIdx() for nbr in mol.GetAtomWithIdx(boron_idx).GetNeighbors()
        if not nbr.HasProp("bsubpc_label")
    )
    bond = mol.GetBondBetweenAtoms(boron_idx, axial_idx)
    return FragmentOnBonds(mol, [bond.GetIdx()], addDummies=True)


def draw_split_bsubpc(drawer, mol: Mol, legend: str = "") -> None:
    """Render a labeled BsubPc as bowl + axial fragment into the given
    MolDraw2D drawer. Caller creates and configures the drawer (size,
    legend font, etc.); after this returns the drawer is finished and
    its output can be retrieved via GetDrawingText()."""
    # AddHs first so that an axial-H BsubPc (no explicit H on boron in
    # the input) still has an atom for the splitter to find.
    mol = AddHs(mol)
    fragged = split_subpc_at_apical_bond(mol)
    draw_mol = RemoveHs(fragged)
    rdDepictor.Compute2DCoords(draw_mol)
    prepared = rdMolDraw2D.PrepareMolForDrawing(draw_mol, kekulize=False)

    options = drawer.drawOptions()
    options.clearBackground = False
    # FragmentOnBonds(addDummies=True) marks the cleavage with isotope-
    # labeled `[*]` atoms; these flags style them as plain attachment
    # stubs instead of map-numbered dummies.
    options.dummiesAreAttachments = True
    options.dummyIsotopeLabels = False

    drawer.DrawMolecule(prepared, legend=legend)
    drawer.FinishDrawing()


def write_diagram(mol: Mol, outpath, width: int = 450, height: int = 320) -> None:
    """Write a single SVG of a labeled BsubPc, split at the boron-axial bond."""
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    draw_split_bsubpc(drawer, mol)
    Path(outpath).write_text(drawer.GetDrawingText(), encoding="utf-8")
