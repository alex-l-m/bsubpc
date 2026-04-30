'''Search the Cambridge Crystallography Database for BsubPc's and save them as mol2 and CIF files.'''
import math
from pathlib import Path
import pandas as pd
import ccdc.io
import ccdc.search

# Template of a BsubPc
template_smarts_raw = '[#5]12-[#7]3:[#6]4:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6]:3-[#7]=[#6]3:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6](:[#7]:3-1)=[#7]-[#6]1=[#7]~2-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#7]-4'
# Make it less specific about bond types
template_smarts = template_smarts_raw.replace(':', '~').replace('=', '~').replace('-', '~')
template_substructure = ccdc.search.SMARTSSubstructure(template_smarts)
MIN_CRYSTAL_DISTANCE = 0.5

# Create a search based on the smarts string
subpc_search = ccdc.search.SubstructureSearch()
subpc_search.add_substructure(template_substructure)

# Execute search and save mol2 files, CIF files, and atom mappings.
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
    """Return molecule atoms after checks required for writing a reliable mapping."""
    atoms = list(molecule.atoms)
    if not molecule.all_atoms_have_sites:
        siteless = [atom.label for atom in atoms if atom.coordinates is None]
        raise ValueError(f'molecule contains atoms without crystal sites: {", ".join(siteless[:10])}')

    molecule_labels = [atom.label for atom in atoms]
    if len(set(molecule_labels)) != len(molecule_labels):
        raise ValueError('molecule atom labels are not unique')
    return atoms


def _cif_atom_site_labels(cif_path: Path) -> list[str]:
    """Return CIF atom-site labels using CCDC's structured CIF reader."""
    entries = ccdc.io.EntryReader(str(cif_path))
    if len(entries) != 1:
        raise ValueError(f'expected one CIF data block in {cif_path}, found {len(entries)}')
    labels = entries[0].attributes['_atom_site_label']
    if not labels:
        raise ValueError(f'{cif_path} does not contain _atom_site_label values')
    return list(labels)


def _molecule_to_cif_mapping(atoms, cif_path: Path) -> list[dict]:
    """Build molecule-index to CIF atom-site-index rows using CCDC data objects."""
    cif_labels = _cif_atom_site_labels(cif_path)
    if len(set(cif_labels)) != len(cif_labels):
        raise ValueError('CIF atom-site labels are not unique')

    label_to_cif_idx = {label: index for index, label in enumerate(cif_labels)}
    missing = [atom.label for atom in atoms if atom.label not in label_to_cif_idx]
    if missing:
        raise ValueError(f'CIF is missing molecule atom labels: {", ".join(missing[:10])}')

    return [
        {
            'molecule_atom_idx': atom_index,
            'crystal_atom_idx': label_to_cif_idx[atom.label],
            'atom_label': atom.label,
            'element': atom.atomic_symbol,
        }
        for atom_index, atom in enumerate(atoms)
    ]


def _base_atom_site_label(label: str) -> str:
    """Strip CCDC's numeric P1 expansion suffix from an atom-site label."""
    base, separator, suffix = label.rpartition('_')
    if separator and suffix.isdigit():
        return base
    return label


def _validate_cif_matches_molecule(molecule_atoms, cif_path: Path) -> None:
    """Require the written CIF atom sites to contain the selected crystal molecule."""
    molecule_labels = [atom.label for atom in molecule_atoms]
    molecule_label_set = set(molecule_labels)
    cif_labels = _cif_atom_site_labels(cif_path)
    cif_base_labels = {_base_atom_site_label(label) for label in cif_labels}
    if not molecule_label_set.issubset(cif_base_labels):
        missing_from_cif = sorted(molecule_label_set - cif_base_labels)
        raise ValueError(
            'written CIF atom sites do not match selected disorder molecule; '
            f'missing from CIF: {", ".join(missing_from_cif[:10])}'
        )
    extra_base_labels = sorted(cif_base_labels - molecule_label_set)
    if extra_base_labels:
        raise ValueError(
            'written CIF atom sites include labels outside the selected disorder molecule; '
            f'extra labels: {", ".join(extra_base_labels[:10])}'
        )


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


_clear_previous_outputs()
seen = set()
outrows = []
mapping_rows = []
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
        if template_substructure.nmatch_molecule(molecule) != 1:
            raise ValueError('matched component does not contain exactly one BsubPc template match')

        atoms = _validated_atoms(molecule)
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
        _validate_cif_matches_molecule(crystal_atoms, cif_output_path)
        mol_mapping_rows = _molecule_to_cif_mapping(atoms, cif_output_path)
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
    for row in mol_mapping_rows:
        mapping_rows.append({'mol_id': mol_id, **row})

outdf = pd.DataFrame(outrows)
outdf.to_csv('search_results.csv', index=False)
pd.DataFrame(
    mapping_rows,
    columns=['mol_id', 'molecule_atom_idx', 'crystal_atom_idx', 'atom_label', 'element'],
).to_csv('atom_mapping.csv', index=False)
pd.DataFrame(exclusion_rows, columns=['mol_id', 'reason']).to_csv('search_exclusions.csv', index=False)
