'''Search the Cambridge Crystallography Database for BsubPc's and save them as mol2 files.'''
from pathlib import Path
import pandas as pd
import ccdc.io
import ccdc.search

# Template of a BsubPc
template_smarts_raw = '[#5]12-[#7]3:[#6]4:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6]:3-[#7]=[#6]3:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6](:[#7]:3-1)=[#7]-[#6]1=[#7]~2-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#7]-4'
# Make it less specific about bond types
template_smarts = template_smarts_raw.replace(':', '~').replace('=', '~').replace('-', '~')
template_substructure = ccdc.search.SMARTSSubstructure(template_smarts)

# Create a search based on the smarts string
subpc_search = ccdc.search.SubstructureSearch()
subpc_search.add_substructure(template_substructure)

# Execute search and save mol2 files and cif files
mol2_output_dir = Path('csd_molecules')
mol2_output_dir.mkdir(exist_ok=True)
cif_output_dir = Path('csd_crystals')
cif_output_dir.mkdir(exist_ok=True)
hits = subpc_search.search()

seen = set()
outrows = []
for hit in hits:
    mol_id = hit.identifier
    if mol_id in seen:
        continue
    seen.add(mol_id)
    mol2_output_path = mol2_output_dir / f'{mol_id}.mol2'
    cif_output_path = cif_output_dir / f'{mol_id}.cif'
    molecule = hit.match_components()[0]
    if template_substructure.nmatch_molecule(molecule) != 1:
            continue
    # Save as mol2
    ccdc.io.MoleculeWriter(mol2_output_path).write(molecule)
    # Save as cif
    ccdc.io.CrystalWriter(cif_output_path).write(hit.crystal)

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

outdf = pd.DataFrame(outrows)
outdf.to_csv('search_results.csv', index=False)
