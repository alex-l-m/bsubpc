'''Search the Cambridge Crystallography Database for BsubPc's and save them as mol2 files'''
from os import mkdir
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

# Execute search and save mol2 files
mol2_output_dir_name = 'search_results'
try:
    mkdir(mol2_output_dir_name)
except FileExistsError:
    pass
hits = subpc_search.search()
print(f'Found {len(hits)} hits')

seen = set()
outrows = []
for hit in hits:
    mol_id = hit.identifier
    if mol_id in seen:
        continue
    seen.add(mol_id)
    mol2_output_path = f'{mol2_output_dir_name}/{mol_id}.mol2'
    molecule = hit.match_components()[0]
    if template_substructure.nmatch_molecule(molecule) != 1:
            print(f"Skipping {mol_id}: != 1 match in extracted component")
            continue
    # Save as mol2
    ccdc.io.MoleculeWriter(mol2_output_path).write(molecule)
    print(f'Wrote {mol_id} to {mol2_output_path}')

    formal_charge = molecule.formal_charge

    outrow = {'mol_id': mol_id, 'formal_charge': formal_charge, 'mol2_path': mol2_output_path}
    outrows.append(outrow)

pd.DataFrame(outrows).to_csv('search_results.csv', index=False)
