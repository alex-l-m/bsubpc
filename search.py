'''Search the Cambridge Crystallography Database for BsubPc's and save them as mol2 and CIF files.'''
import math
from pathlib import Path
import pandas as pd
import ccdc.io
import ccdc.search

from bsubpc_match import TEMPLATE_ATOM_LABELS, TEMPLATE_SMARTS


template_substructure = ccdc.search.SMARTSSubstructure(TEMPLATE_SMARTS)
MIN_CRYSTAL_DISTANCE = 0.5
MOL2_MATCH_PATH = Path('bsubpc_mol2_match_indices.csv')
CIF_MATCH_PATH = Path('bsubpc_cif_match_indices.csv')

# Create a search based on the smarts string
subpc_search = ccdc.search.SubstructureSearch()
subpc_search.add_substructure(template_substructure)

# Execute search and save mol2 files, CIF files, and match index tables.
mol2_output_dir = Path('csd_molecules')
mol2_output_dir.mkdir(exist_ok=True)
cif_output_dir = Path('csd_crystals')
cif_output_dir.mkdir(exist_ok=True)
hits = subpc_search.search()


def _clear_previous_outputs() -> None:
    """Remove stale generated outputs before writing a new filtered result set."""
    paths = [
        *mol2_output_dir.glob('*.mol2'),
        *cif_output_dir.glob('*.cif'),
    ]
    for path in paths:
        path.unlink()


def _validated_atoms(molecule):
    """Return molecule atoms after checks required for writing reliable outputs."""
    atoms = list(molecule.atoms)
    if not molecule.all_atoms_have_sites:
        siteless = [atom.label for atom in atoms if atom.coordinates is None]
        raise ValueError(f'molecule contains atoms without crystal sites: {", ".join(siteless[:10])}')

    molecule_labels = [atom.label for atom in atoms]
    if len(set(molecule_labels)) != len(molecule_labels):
        raise ValueError('molecule atom labels are not unique')
    return atoms


def _validate_no_disordered_heavy_atoms(crystal) -> None:
    """Reject crystal structures with any non-hydrogen disorder-group atoms."""
    if not crystal.has_disorder:
        return

    disordered_heavy_atom_labels = set()
    for assembly in crystal.disorder.assemblies:
        for group in assembly.groups:
            disordered_heavy_atom_labels.update(
                atom.label
                for atom in group.atoms
                if atom.atomic_symbol != 'H'
            )
    if disordered_heavy_atom_labels:
        labels = ', '.join(sorted(disordered_heavy_atom_labels)[:10])
        raise ValueError(f'crystal has disordered heavy atoms: {labels}')


def _validate_min_distance(atoms) -> None:
    """Reject active crystal contents with impossible close contacts."""
    min_distance = math.inf
    closest_labels = None
    for i, atom_i in enumerate(atoms):
        pos_i = atom_i.coordinates
        for atom_j in atoms[i + 1:]:
            pos_j = atom_j.coordinates
            distance = math.dist((pos_i.x, pos_i.y, pos_i.z), (pos_j.x, pos_j.y, pos_j.z))
            if distance < min_distance:
                min_distance = distance
                closest_labels = (atom_i.label, atom_j.label)
    if min_distance < MIN_CRYSTAL_DISTANCE:
        raise ValueError(
            f'active crystal atoms {closest_labels[0]} and {closest_labels[1]} '
            f'are only {min_distance:.3f} A apart'
        )


def _bsubpc_label(match_index: int) -> str:
    """Return the optional human-readable label used by match.py."""
    return TEMPLATE_ATOM_LABELS.get(match_index, '')


def _single_component_match_atoms(molecule):
    """Return the unique BsubPc match in a molecule using CCDC's matcher."""
    # Use full substructure hits here, not template_substructure.nmatch_molecule().
    # nmatch_molecule() counts atoms that can match the query start atom; it is
    # not the number of complete BsubPc substructure matches.
    matches = subpc_search.search(molecule)
    if len(matches) != 1:
        raise ValueError(f'matched component contains {len(matches)} BsubPc template matches')
    # CCDC's SubstructureHit.match_atoms() follows the motif atom_match order,
    # which is the same query-atom order used by RDKit's GetSubstructMatches().
    # Enumerating this list therefore gives the bsubpc_idx values used by match.py.
    return list(matches[0].match_atoms())


def _mol2_match_rows(mol_id: str, match_atoms) -> list[dict]:
    """Build match-index rows for the written mol2 component."""
    return [
        {
            'mol_id': mol_id,
            'bsubpc_idx': match_index,
            'bsubpc_label': _bsubpc_label(match_index),
            'mol2_atom_idx': atom.index,
            'atom_label': atom.label,
            'element': atom.atomic_symbol,
        }
        for match_index, atom in enumerate(match_atoms)
    ]


def _cif_match_rows(mol_id: str, cif_path: Path) -> list[dict]:
    """Build match-index rows for every BsubPc copy in the written CIF."""
    # Search the written CIF path, rather than the in-memory crystal, so CCDC
    # reads exactly the atom-site order and unique labels emitted by
    # CrystalWriter. This matters after P1 reduction: repeated atom labels in
    # the crystal object are written as unique labels such as B1_2, B1_3, etc.
    matches = [list(hit.match_atoms()) for hit in subpc_search.search(str(cif_path))]
    if not matches:
        raise ValueError('written CIF contains no BsubPc template matches')
    # The P1 CIF can contain several BsubPc copies. Sort by the matched boron
    # atom index so cif_match_idx is stable and follows the CIF atom-site order.
    matches.sort(key=lambda atoms: atoms[0].index)

    rows = []
    for cif_match_idx, match_atoms in enumerate(matches):
        for match_index, atom in enumerate(match_atoms):
            rows.append(
                {
                    'mol_id': mol_id,
                    'cif_match_idx': cif_match_idx,
                    'bsubpc_idx': match_index,
                    'bsubpc_label': _bsubpc_label(match_index),
                    'cif_atom_idx': atom.index,
                    'atom_site_label': atom.label,
                    'element': atom.atomic_symbol,
                }
            )
    return rows


_clear_previous_outputs()
seen = set()
outrows = []
mol2_match_rows = []
cif_match_rows = []
exclusion_rows = []
for hit in hits:
    mol_id = hit.identifier
    if mol_id in seen:
        continue
    seen.add(mol_id)
    mol2_output_path = mol2_output_dir / f'{mol_id}.mol2'
    cif_output_path = cif_output_dir / f'{mol_id}.cif'
    try:
        molecule = hit.match_components()[0]
        # Rematch the component that will be written to mol2 so mol2_atom_idx
        # refers to the output molecule's atom order, not the parent CSD entry.
        component_match_atoms = _single_component_match_atoms(molecule)
        _validated_atoms(molecule)
        crystal = hit.crystal
        _validate_no_disordered_heavy_atoms(crystal)
        # crystal.molecule is the CSD-selected disorder state, not necessarily
        # the same as the crystal's editable molecule used by CrystalWriter.
        crystal_molecule = crystal.molecule
        crystal_atoms = _validated_atoms(crystal_molecule)
        _validate_min_distance(crystal_atoms)

        # Save as mol2.
        ccdc.io.MoleculeWriter(mol2_output_path).write(molecule)
        # This setter is counterintuitive but necessary: it materializes the
        # selected disorder molecule into the editable crystal structure.
        # Without it, CrystalWriter can emit the full disordered atom-site table.
        # Reduce to P1 after that selection so the CIF already contains the
        # same full unit-cell expansion that downstream extxyz generation uses.
        crystal.molecule = crystal_molecule
        crystal = crystal.reduce_symmetry_to_p1()
        ccdc.io.CrystalWriter(cif_output_path).write(crystal)
        mol2_match_rows_for_hit = _mol2_match_rows(mol_id, component_match_atoms)
        cif_match_rows_for_hit = _cif_match_rows(mol_id, cif_output_path)
    except Exception as exc:
        if mol2_output_path.exists():
            mol2_output_path.unlink()
        if cif_output_path.exists():
            cif_output_path.unlink()
        exclusion_rows.append(
            {
                'mol_id': mol_id,
                'reason': f'{type(exc).__name__}: {exc}',
            }
        )
        continue

    formal_charge = molecule.formal_charge

    outrow = {
        # Original search_results.csv columns.
        'mol_id': mol_id,
        'formal_charge': formal_charge,
        'mol2_path': str(mol2_output_path),
        'cif_path': str(cif_output_path),
        # Compatibility columns for oled-database scripts.  In this search, the
        # CSD hit identifier is also the structure entry/refcode that the
        # extractor should read back from the CSD.
        'structure_entry': mol_id,
        'structure_path': str(mol2_output_path),
        'csd_refcode': mol_id,
        'entry': mol_id,
    }
    outrows.append(outrow)
    mol2_match_rows.extend(mol2_match_rows_for_hit)
    cif_match_rows.extend(cif_match_rows_for_hit)

outdf = pd.DataFrame(outrows)
outdf.to_csv('search_results.csv', index=False)
pd.DataFrame(
    mol2_match_rows,
    columns=['mol_id', 'bsubpc_idx', 'bsubpc_label', 'mol2_atom_idx', 'atom_label', 'element'],
).to_csv(MOL2_MATCH_PATH, index=False)
pd.DataFrame(
    cif_match_rows,
    columns=[
        'mol_id',
        'cif_match_idx',
        'bsubpc_idx',
        'bsubpc_label',
        'cif_atom_idx',
        'atom_site_label',
        'element',
    ],
).to_csv(CIF_MATCH_PATH, index=False)
pd.DataFrame(exclusion_rows, columns=['mol_id', 'reason']).to_csv('search_exclusions.csv', index=False)
