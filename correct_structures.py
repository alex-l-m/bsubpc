from glob import glob
from os.path import basename, splitext
from pathlib import Path
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdMolTransforms import CanonicalizeConformer
from rdkit.Chem.rdmolfiles import MolFromMol2File, MolToMolFile, MolFromSmarts

indir = 'search_results'
inpaths = glob(f'{indir}/*.mol2')

outdir = Path('crystal_geom')
outdir.mkdir(exist_ok=True)

# Template of a BsubPc
template_smarts_raw = '[#5]12-[#7]3:[#6]4:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6]:3-[#7]=[#6]3:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6](:[#7]:3-1)=[#7]-[#6]1=[#7]~2-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#7]-4'
# Make it less specific about bond types
template_smarts = template_smarts_raw.replace(':', '~').replace('=', '~').replace('-', '~')
template = MolFromSmarts(template_smarts)

def single_match(mol: Mol) -> bool:
    '''Return true if there's only one match to the substructure. Return false
    if the substructure can be matched to non-overlapping sets of atoms.'''
    matches = mol.GetSubstructMatches(template)
    first_match = matches[0]
    for match in matches:
        if set(match) != set(first_match):
            return False
    return True

def translate_rotate(mol: Mol) -> None:
    '''Move the molecule to the origin and rotate. Currently just uses
    CanonicalizeConformer, I might make a custom version later. Assumes there's
    only one conformer'''
    conformer = mol.GetConformer(0)
    CanonicalizeConformer(conformer)

for inpath in inpaths:
    mol = MolFromMol2File(inpath)
    if mol is None:
        print(f'Error reading {inpath}')
        continue
    if not mol.HasSubstructMatch(template):
        print(f'No match to template in {inpath}')
        continue
    if not single_match(mol):
        print(f'Multiple matches to template in {inpath}')
        continue
    mol_id, ext = splitext(basename(inpath))
    mol.SetProp('_Name', mol_id)
    translate_rotate(mol)
    outpath = outdir / f'{mol_id}.mol'
    MolToMolFile(mol, outpath)
