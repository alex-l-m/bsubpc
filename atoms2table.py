'''Given a ASE atoms object as input, output a table containing both positions and forces'''
import sys
import os.path
import pandas as pd
import ase.io

inpath = sys.argv[1]
stem, ext = os.path.splitext(inpath)
outpath = f'{stem}.csv'

first = True
for i, atoms in enumerate(ase.io.iread(inpath, index=':', format='extxyz')):
    positions = atoms.get_positions()
    forces = atoms.get_forces()
    df = pd.DataFrame({
        'frame': [i] * len(atoms),
        'atom_id': range(len(atoms)),
        'element': atoms.get_chemical_symbols(),
        'x': positions[:, 0],
        'y': positions[:, 1],
        'z': positions[:, 2],
        'fx': forces[:, 0],
        'fy': forces[:, 1],
        'fz': forces[:, 2]
    })
    df.to_csv(outpath, index=False, mode='w' if first else 'a', header=first)
    first = False
