from pathlib import Path
from rdkit.Chem import rdDepictor
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolops import FragmentOnBonds, RemoveHs, AddHs

from match import assert_labeled

def _split_subpc_at_apical_bond(mol: Mol) -> Mol:
    """Return (axial_substituent_fragment, remainder_fragment)."""
    assert_labeled(mol)
    for atom in mol.GetAtoms():
        if atom.HasProp("bsubpc_label") and \
                atom.GetProp("bsubpc_label") == "boron":
            boron_idx = atom.GetIdx()
            break
    # Find the unlabeled neighbor
    for nbr in mol.GetAtomWithIdx(boron_idx).GetNeighbors():
        if not nbr.HasProp("bsubpc_label"):
            axial_neighbor_idx = nbr.GetIdx()
            break
    # Fragment on the bond
    broken_bond = mol.GetBondBetweenAtoms(boron_idx, axial_neighbor_idx)
    fragged = FragmentOnBonds(mol, [broken_bond.GetIdx()], addDummies=True)
    return fragged

def _draw_svg(mol: Mol, out_path: Path, width: int = 450, height: int = 320) -> None:
    draw_mol = RemoveHs(mol)
    rdDepictor.Compute2DCoords(draw_mol)
    prepared = rdMolDraw2D.PrepareMolForDrawing(draw_mol, kekulize=False)

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    options.addStereoAnnotation = False
    options.clearBackground = False
    options.dummiesAreAttachments = True
    options.dummyIsotopeLabels = False

    drawer.DrawMolecule(prepared)
    drawer.FinishDrawing()
    out_path.write_text(drawer.GetDrawingText(), encoding="utf-8")

def write_diagram(mol: Mol, outpath) -> None:
    '''Write SVG of axial substituent and BsubPc bowl to separate files'''
    # Need to add hydrogens or there will be an issue when the axial
    # substituent is just hydrogen
    mol = AddHs(mol)
    fragged = _split_subpc_at_apical_bond(mol)
    _draw_svg(fragged, Path(outpath))
