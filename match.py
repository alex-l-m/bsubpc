"""Match and label the BsubPc core atoms."""
import re
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdmolfiles import MolFromSmarts


TEMPLATE_SMARTS = (
    "[#5]12-[#7]3:[#6]4:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6]:3-"
    "[#7]=[#6]3:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6](:[#7]:3-1)="
    "[#7]-[#6]1=[#7]~2-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#7]-4"
).replace(":", "~").replace("=", "~").replace("-", "~")
TEMPLATE = MolFromSmarts(TEMPLATE_SMARTS)


TEMPLATE_ATOM_LABELS = {
    0: "boron",
    1: "pyrrole_nitrogen_1",
    19: "pyrrole_nitrogen_2",
    22: "pyrrole_nitrogen_3",
    10: "imine_nitrogen_1",
    20: "imine_nitrogen_2",
    30: "imine_nitrogen_3",
    6: "outer_terminal_carbon_1",
    7: "outer_terminal_carbon_2",
    15: "outer_terminal_carbon_3",
    16: "outer_terminal_carbon_4",
    27: "outer_terminal_carbon_5",
    28: "outer_terminal_carbon_6",
}


CATEGORY_REGEX = re.compile(r"^(.*)_\d+$")


TEMPLATE_ATOM_INDEX = {label: index for index, label in TEMPLATE_ATOM_LABELS.items()}


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
