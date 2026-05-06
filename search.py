'''Search the Cambridge Crystallography Database for BsubPc's and save them as mol2 and CIF files.

For each CSD hit, the matched-component mol2 and the P1-expanded crystal CIF
are written **independently**: a problem that affects only the crystal (e.g.
unphysically close atoms in the unit cell) no longer poisons the per-molecule
mol2 that downstream sTDA work consumes, and vice versa.

`search_results.csv` has one row per hit that produced **either** output. Empty
`mol2_path` / `cif_path` cells indicate the corresponding output was rejected;
`mol2_status` / `cif_status` give per-output ok/skipped/failed and their
reasons live in `mol2_error` / `cif_error`. The matched component's element-
count signature ends up in `mol2_formula`, so `process_search_results.py` can
gate the mol output against the entry's CSD chemical diagram without
re-parsing the mol2. `search_exclusions.csv` keeps a row-per-failure log for
quick inspection.
'''
import argparse
import math
from collections import Counter
from pathlib import Path
import pandas as pd
import ccdc.io
import ccdc.search

from bsubpc_match import TEMPLATE_ATOM_LABELS, TEMPLATE_SMARTS


template_substructure = ccdc.search.SMARTSSubstructure(TEMPLATE_SMARTS)
MIN_CRYSTAL_DISTANCE = 0.5
MOL2_MATCH_PATH = Path('bsubpc_mol2_match_indices.csv')
CIF_MATCH_PATH = Path('bsubpc_cif_match_indices.csv')


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--prefix-filter',
        help='Only search CSD refcodes starting with this prefix; intended for quick testing.',
    )
    return parser.parse_args()


def _new_subpc_search():
    """Create the BsubPc substructure search used for CSD and local matching."""
    search = ccdc.search.SubstructureSearch()
    search.add_substructure(template_substructure)
    # CCDC's no_disorder = 'Non-hydrogen' database filter rejects entries that
    # have any disorder on heavy atoms — even when the active disorder state
    # selected by CCDC is fully populated and the matched BsubPc component is
    # clean. We've seen it drop refcodes such as XOGYEZ and ELIQUO that we
    # absolutely want in the sTDA set, so the filter stays off here. Per-hit
    # validation in `_validated_atoms` still rejects components with
    # unresolved (siteless) atoms.
    return search


subpc_search = _new_subpc_search()


def _search_hits(prefix_filter: str | None):
    """Run the CSD search, optionally constraining the CSD query by refcode prefix."""
    search = _new_subpc_search()
    if prefix_filter:
        prefix_search = ccdc.search.TextNumericSearch()
        # TextNumericSearch supports STARTS_WITH via mode='start', so the
        # prefix restriction is applied in the CSD query rather than after the
        # expensive substructure search.
        prefix_search.add_identifier(prefix_filter.upper(), mode='start')
        search = ccdc.search.CombinedSearch(search & prefix_search)
    return search.search()

# Execute search and save mol2 files, CIF files, and match index tables.
mol2_output_dir = Path('csd_molecules')
mol2_output_dir.mkdir(exist_ok=True)
cif_output_dir = Path('csd_crystals')
cif_output_dir.mkdir(exist_ok=True)


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


def _validate_min_distance(atoms) -> None:
    """Reject active crystal contents with impossible close contacts."""
    # This is an absolute coordinate-distance cutoff. Crystal.contacts() is a
    # VdW-corrected nonbonded contact search, so it cannot express this exact
    # sub-angstrom sanity check directly.
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


def _format_error(exc: BaseException) -> str:
    return f'{type(exc).__name__}: {exc}'


def _formula_string(molecule) -> str:
    """Canonical sorted element-count signature for a CCDC molecule."""
    counter = Counter(a.atomic_symbol for a in molecule.atoms)
    return ''.join(f'{el}{counter[el]}' for el in sorted(counter))


def _process_mol2(mol_id: str, hit, mol2_output_path: Path):
    """Try to write the matched-component mol2.

    Returns ``(status, error, match_rows, molecule, formula)``. ``formula`` is
    the matched component's sorted element-count signature; downstream uses it
    to verify the mol2 against the CSD chemical-diagram view formulas.
    """
    try:
        molecule = hit.match_components()[0]
        component_match_atoms = _single_component_match_atoms(molecule)
        _validated_atoms(molecule)
        ccdc.io.MoleculeWriter(mol2_output_path).write(molecule)
    except Exception as exc:
        if mol2_output_path.exists():
            mol2_output_path.unlink()
        return 'failed', _format_error(exc), [], None, ''
    return (
        'ok',
        '',
        _mol2_match_rows(mol_id, component_match_atoms),
        molecule,
        _formula_string(molecule),
    )


def _process_cif(mol_id: str, hit, cif_output_path: Path):
    """Try to write the P1-expanded crystal cif. Returns (status, error, match_rows)."""
    try:
        crystal = hit.crystal
        # crystal.molecule is the CSD-selected disorder state, not necessarily
        # the same as the crystal's editable molecule used by CrystalWriter.
        crystal_molecule = crystal.molecule
        crystal_atoms = _validated_atoms(crystal_molecule)
        _validate_min_distance(crystal_atoms)
        # This setter is counterintuitive but necessary: it materializes the
        # selected disorder molecule into the editable crystal structure.
        # Without it, CrystalWriter can emit the full disordered atom-site table.
        # Reduce to P1 after that selection so the CIF already contains the
        # same full unit-cell expansion that downstream extxyz generation uses.
        crystal.molecule = crystal_molecule
        crystal = crystal.reduce_symmetry_to_p1()
        ccdc.io.CrystalWriter(cif_output_path).write(crystal)
        rows = _cif_match_rows(mol_id, cif_output_path)
    except Exception as exc:
        if cif_output_path.exists():
            cif_output_path.unlink()
        return 'failed', _format_error(exc), []
    return 'ok', '', rows


args = _parse_args()
hits = _search_hits(args.prefix_filter)
_clear_previous_outputs()
seen_deposition_numbers = set()
outrows = []
mol2_match_rows = []
cif_match_rows = []
exclusion_rows = []
for hit in hits:
    mol_id = hit.identifier
    mol2_output_path = mol2_output_dir / f'{mol_id}.mol2'
    cif_output_path = cif_output_dir / f'{mol_id}.cif'

    entry = hit.entry
    deposition_number = entry.ccdc_number
    # Deposition number is distinct from the six-letter CSD refcode, and is
    # sometimes associated with multiple entries with different refcodes.
    # De-duplicate after at least one successful output, not before, so a
    # failed refcode does not suppress a later related hit that would work.
    if deposition_number is not None and deposition_number in seen_deposition_numbers:
        continue

    mol2_status, mol2_error, mol2_rows_for_hit, molecule, mol2_formula = _process_mol2(
        mol_id, hit, mol2_output_path
    )
    cif_status, cif_error, cif_rows_for_hit = _process_cif(
        mol_id, hit, cif_output_path
    )

    if mol2_status != 'ok':
        exclusion_rows.append({'mol_id': mol_id, 'output': 'mol2', 'reason': mol2_error})
    if cif_status != 'ok':
        exclusion_rows.append({'mol_id': mol_id, 'output': 'cif', 'reason': cif_error})
    if mol2_status != 'ok' and cif_status != 'ok':
        # Hit produced no usable output; skip it from search_results.csv too.
        continue

    # Determine formal charge from whichever output we got. Prefer the
    # matched-component charge since that is what downstream mol2 readers see.
    if molecule is not None:
        formal_charge = molecule.formal_charge
    else:
        formal_charge = hit.crystal.molecule.formal_charge

    outrows.append(
        {
            'mol_id': mol_id,
            'ccdc_number': deposition_number,
            'formal_charge': formal_charge,
            'mol2_path': str(mol2_output_path) if mol2_status == 'ok' else '',
            'cif_path': str(cif_output_path) if cif_status == 'ok' else '',
            'mol2_status': mol2_status,
            'cif_status': cif_status,
            'mol2_error': mol2_error,
            'cif_error': cif_error,
            'mol2_formula': mol2_formula,
            # Compatibility columns for oled-database scripts.  In this search,
            # the CSD hit identifier is also the structure entry/refcode that
            # the extractor should read back from the CSD.
            'structure_entry': mol_id,
            'structure_path': str(mol2_output_path) if mol2_status == 'ok' else '',
            'csd_refcode': mol_id,
            'entry': mol_id,
        }
    )
    mol2_match_rows.extend(mol2_rows_for_hit)
    cif_match_rows.extend(cif_rows_for_hit)
    if deposition_number is not None:
        seen_deposition_numbers.add(deposition_number)

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
pd.DataFrame(exclusion_rows, columns=['mol_id', 'output', 'reason']).to_csv(
    'search_exclusions.csv', index=False
)
