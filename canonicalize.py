"""Translate and rotate a BsubPc into a canonical orientation.

The canonical frame is defined by three geometric constraints:

* the boron sits at the origin;
* the bowl lies below the xy plane;
* the plane through boron, +z, and ``pyrrole_nitrogen_1`` is the xz plane.

The vertical axis is the normal to the best-fit plane of the matched BsubPc core,
obtained from the smallest right-singular vector of the centered core coordinates.
No explicit symmetry detection is used.
"""

import numpy as np
from rdkit.Chem.rdchem import Mol
from scipy.spatial.transform import Rotation

from match import label_core_atoms, assert_labeled


VERTICAL_AXIS = np.array([0.0, 0.0, 1.0])
REFERENCE_LABEL = "pyrrole_nitrogen_1"
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
    """Return the unit normal to the best-fit plane of centered positions."""
    _, _, vh = np.linalg.svd(np.asarray(pos, dtype=float), full_matrices=False)
    return vh[-1] / np.linalg.norm(vh[-1])


def rotate_to_vertical(mol: Mol) -> None:
    """Rotate every conformer so that the BsubPc bowl lies below the xy plane."""
    assert_centered(mol)
    core_idx = _core_atom_indices(mol)

    for conformer in mol.GetConformers():
        pos = conformer.GetPositions()
        axis = vertical_axis_from_pos(pos[core_idx])
        if np.mean(pos[core_idx] @ axis) > 0:
            axis = -axis
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


def canonicalize_bsubpc(mol: Mol) -> None:
    """Label, center, and rotate ``mol`` into the canonical BsubPc frame in-place."""
    label_core_atoms(mol)
    center_on_boron(mol)
    rotate_to_vertical(mol)
    rotate_around_vertical(mol)
