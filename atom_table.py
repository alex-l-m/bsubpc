from typing import Iterable
import pandas as pd
from rdkit.Chem.rdchem import Mol
from match import assert_labeled

def atom_frame(mol: Mol) -> pd.DataFrame:
    assert_labeled(mol)
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
            "bsubpc_idx": [atom.GetProp("bsubpc_idx") if atom.HasProp("bsubpc_idx") else None for atom in atoms],
            "bsubpc_label": [atom.GetProp("bsubpc_label") if atom.HasProp("bsubpc_label") else None for atom in atoms],
        }
    )

def mols2tbl(mols: Iterable[Mol]) -> pd.DataFrame:
    return pd.concat(map(atom_frame, mols), ignore_index=True)
