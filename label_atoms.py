"""Label the BsubPc core atoms used for bowl-depth analysis and export coordinates."""

from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.Draw import MolToImage

RDLogger.DisableLog("rdApp.warning")

MOL_DIR = Path("crystal_geom")
OUTPUT_PREFIX = "crystal"

TEMPLATE_SMARTS = (
    "[#5]12-[#7]3:[#6]4:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6]:3-"
    "[#7]=[#6]3:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6](:[#7]:3-1)="
    "[#7]-[#6]1=[#7]~2-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#7]-4"
).replace(":", "~").replace("=", "~").replace("-", "~")
TEMPLATE = Chem.MolFromSmarts(TEMPLATE_SMARTS)

TEMPLATE_INDEX_CATEGORY = {
    0: "boron",
    1: "pyrrole_nitrogen",
    19: "pyrrole_nitrogen",
    22: "pyrrole_nitrogen",
    10: "imine_nitrogen",
    20: "imine_nitrogen",
    30: "imine_nitrogen",
    6: "outer_terminal_carbon",
    7: "outer_terminal_carbon",
    15: "outer_terminal_carbon",
    16: "outer_terminal_carbon",
    27: "outer_terminal_carbon",
    28: "outer_terminal_carbon",
}


def load_mol(path: Path):
    readers = {
        ".mol": lambda p: Chem.MolFromMolFile(str(p), removeHs=True),
        ".mol2": lambda p: Chem.MolFromMol2File(str(p), removeHs=True),
    }
    mol = readers[path.suffix](path)
    if mol is None:
        raise ValueError(f"Could not read {path}")
    mol.SetProp("_Name", path.stem)
    return mol


def label_core_atoms(mol):
    match = mol.GetSubstructMatch(TEMPLATE)
    if not match:
        raise ValueError(f"{mol.GetProp('_Name')} does not match the BsubPc template")

    category_by_atom_index = {match[i]: category for i, category in TEMPLATE_INDEX_CATEGORY.items()}
    for atom in mol.GetAtoms():
        category = category_by_atom_index.get(atom.GetIdx(), "")
        atom.SetProp("category", category)
        atom.SetProp("atomNote", category)
    return mol


def atom_frame(mol):
    atoms = list(mol.GetAtoms())
    xyz = mol.GetConformer().GetPositions()
    return pd.DataFrame(
        {
            "mol_id": mol.GetProp("_Name"),
            "atom_id": [atom.GetIdx() for atom in atoms],
            "atom_name": [atom.GetProp("_TriposAtomName") if atom.HasProp("_TriposAtomName") else None for atom in atoms],
            "symbol": [atom.GetSymbol() for atom in atoms],
            "formal_charge": [atom.GetFormalCharge() for atom in atoms],
            "x": xyz[:, 0],
            "y": xyz[:, 1],
            "z": xyz[:, 2],
            "category": [atom.GetProp("category") for atom in atoms],
        }
    )


def save_template_image(path: Path):
    template = Chem.Mol(TEMPLATE)
    for atom in template.GetAtoms():
        atom.SetProp("molAtomMapNumber", str(atom.GetIdx()))
    MolToImage(template).save(path)


def save_labeled_images(mols, out_dir: Path):
    out_dir.mkdir(exist_ok=True)
    for mol in mols:
        MolToImage(mol).save(out_dir / f"{mol.GetProp('_Name')}.png")


def main():
    mols = [label_core_atoms(load_mol(path)) for path in sorted(MOL_DIR.glob("*.mol*"))]
    save_template_image(Path("template_img.png"))
    save_labeled_images(mols, Path(f"{OUTPUT_PREFIX}_labeled_images"))
    pd.concat(map(atom_frame, mols), ignore_index=True).to_csv(
        f"{OUTPUT_PREFIX}_atom_table.csv", index=False, float_format="%.6f"
    )


if __name__ == "__main__":
    main()
