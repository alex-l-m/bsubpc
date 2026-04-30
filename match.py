"""Match and label the BsubPc core atoms."""
import re
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolFromSmarts

from bsubpc_match import TEMPLATE_ATOM_INDEX, TEMPLATE_ATOM_LABELS, TEMPLATE_SMARTS

TEMPLATE = MolFromSmarts(TEMPLATE_SMARTS)


CATEGORY_REGEX = re.compile(r"^(.*)_\d+$")


def label_core_atoms(mol: Mol) -> None:
    """Annotate the unique BsubPc core match in-place.

    Every matched atom gets a ``bsubpc_idx`` property equal to its query-atom
    index in ``TEMPLATE``. The named atoms in ``TEMPLATE_ATOM_LABELS`` also get a
    human-readable ``bsubpc_label`` property.
    """
    matches = mol.GetSubstructMatches(TEMPLATE)
    if len(matches) != 1:
        raise ValueError(f"expected exactly one BsubPc match, found {len(matches)}")

    for match_index, atom_index in enumerate(matches[0]):
        atom = mol.GetAtomWithIdx(atom_index)
        atom.SetProp("bsubpc_idx", str(match_index))
        if match_index in TEMPLATE_ATOM_LABELS:
            atom.SetProp("bsubpc_label", TEMPLATE_ATOM_LABELS[match_index])


def assert_labeled(mol: Mol) -> None:
    """Assert that ``mol`` carries BsubPc match labels."""
    if not any(atom.HasProp("bsubpc_idx") for atom in mol.GetAtoms()):
        raise ValueError("molecule does not have BsubPc labels; run label_core_atoms first")
