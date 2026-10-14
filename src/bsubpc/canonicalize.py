"""Translate and rotate a BsubPc into a canonical orientation.

The canonical frame is defined by three geometric constraints:

* the boron sits at the origin;
* the bowl lies below the xy plane;
* the plane through boron, +z, and ``pyrrole_nitrogen_1`` is the xz plane;
* ``pyrrole_nitrogen_3`` has positive y.

The vertical axis is the normal to the best-fit plane of the matched BsubPc core,
obtained from the smallest right-singular vector of the centered core coordinates.
The last constraint is met by choosing among the six symmetry-equivalent core
labelings rather than by reflecting, so chirality is preserved.
"""

import numpy as np
from rdkit.Chem.rdchem import Mol
from scipy.spatial.transform import Rotation

from .label import label_core_atoms, assert_labeled
from .strings import TEMPLATE_ATOM_INDEX


VERTICAL_AXIS = np.array([0.0, 0.0, 1.0])
REFERENCE_LABEL = "pyrrole_nitrogen_1"
HANDED_LABEL = "pyrrole_nitrogen_3"
BORON_LABEL = "boron"



def _atom_index_by_label(mol: Mol, label: str) -> int:
    for atom in mol.GetAtoms():
        if atom.HasProp("bsubpc_label") and atom.GetProp("bsubpc_label") == label:
            return atom.GetIdx()
    raise ValueError(f"molecule is missing the labeled atom {label!r}")


def _core_atom_indices(mol: Mol) -> list[int]:
    return [atom.GetIdx() for atom in mol.GetAtoms() if atom.HasProp("bsubpc_idx")]


def _set_positions(conformer, pos: np.ndarray) -> None:
    for atom_index, xyz in enumerate(pos):
        conformer.SetAtomPosition(atom_index, xyz)


def center_on_boron(mol: Mol) -> None:
    """Translate every conformer so that the labeled BsubPc boron is at the origin."""
    assert_labeled(mol)
    boron_idx = _atom_index_by_label(mol, BORON_LABEL)
    for conformer in mol.GetConformers():
        pos = conformer.GetPositions()
        _set_positions(conformer, pos - pos[boron_idx])


def assert_centered(mol: Mol) -> None:
    """Assert that the labeled BsubPc boron sits at the origin in every conformer."""
    assert_labeled(mol)
    boron_idx = _atom_index_by_label(mol, BORON_LABEL)
    for conformer in mol.GetConformers():
        boron_pos = np.asarray(conformer.GetAtomPosition(boron_idx), dtype=float)
        if not np.allclose(boron_pos, np.zeros(3)):
            raise ValueError("molecule is not centered on boron")


def vertical_axis_from_pos(pos: np.ndarray) -> np.ndarray:
    """Return the unit normal to the best-fit plane of centered core positions, pointing away from the bowl."""
    _, _, vh = np.linalg.svd(np.asarray(pos, dtype=float), full_matrices=False)
    axis = vh[-1] / np.linalg.norm(vh[-1])
    return -axis if np.mean(pos @ axis) > 0 else axis


def _is_right_handed(pos: np.ndarray, match: tuple[int, ...]) -> bool:
    """Whether ``match`` labels the core so that N1 -> N3 runs counterclockwise seen from +z."""
    core = pos[list(match)] - pos[match[0]]
    n1 = core[TEMPLATE_ATOM_INDEX[REFERENCE_LABEL]]
    n3 = core[TEMPLATE_ATOM_INDEX[HANDED_LABEL]]
    return np.cross(n1, n3) @ vertical_axis_from_pos(core) > 0


def rotate_to_vertical(mol: Mol) -> None:
    """Rotate every conformer so that the BsubPc bowl lies below the xy plane."""
    assert_centered(mol)
    core_idx = _core_atom_indices(mol)

    for conformer in mol.GetConformers():
        pos = conformer.GetPositions()
        axis = vertical_axis_from_pos(pos[core_idx])
        rotation, _ = Rotation.align_vectors(VERTICAL_AXIS[None, :], axis[None, :])
        _set_positions(conformer, rotation.apply(pos))


def rotate_around_vertical(mol: Mol) -> None:
    """Rotate every conformer about +z so the reference pyrrolic nitrogen lies in the xz plane."""
    assert_centered(mol)
    reference_idx = _atom_index_by_label(mol, REFERENCE_LABEL)

    for conformer in mol.GetConformers():
        pos = conformer.GetPositions()
        reference_xy = np.asarray(pos[reference_idx], dtype=float).copy()
        reference_xy[2] = 0.0
        if np.allclose(reference_xy, 0.0):
            raise ValueError(f"labeled atom {REFERENCE_LABEL!r} lies on the vertical axis")
        angle = np.arctan2(reference_xy[1], reference_xy[0])
        rotation = Rotation.from_rotvec(-angle * VERTICAL_AXIS)
        _set_positions(conformer, rotation.apply(pos))

def canonicalize(mol: Mol) -> None:
    """Label, center, and rotate ``mol`` into the canonical BsubPc frame in-place."""
    pos = mol.GetConformer().GetPositions()
    label_core_atoms(mol, accept=lambda match: _is_right_handed(pos, match))
    center_on_boron(mol)
    rotate_to_vertical(mol)
    rotate_around_vertical(mol)
