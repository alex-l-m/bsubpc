'''Read the STDA log files for each molecules and output a spreadsheet containing the energies'''
import sys
from glob import glob
from os.path import basename, splitext
import pandas as pd
import re

def log2energy(stda_log):
    '''From an stda log as a string, extract and return the excitation energy
    to the lowest excited state as a number, in eV'''
    right_part = False
    for line in stda_log.split("\n"):
        # Based on observation (not inspection of the code) it looks like it
        # always prints three digits after the decimal point. And with
        # sufficiently small energies, there's not necessarily a space between
        # the energy and the wavelength. Therefore the regex can't use * at the
        # end, it has to actually count digits
        match_string = r"\s*1\s*([0-9]*\.[0-9][0-9]?[0-9]?)"
        # Section heading. Newer / large-molecule sTDA runs print the section
        # heading with the singular "excitation energy"; older / small-molecule
        # runs print plural "excitation energies".
        if (
            "excitation energies, transition moments and TDA amplitudes" in line
            or "excitation energy, transition moments and TDA amplitudes" in line
        ):
            right_part = True
            continue
        if line.strip() == "":
            continue
        energy_match = re.match(match_string, line)
        if right_part and energy_match is not None:
            return float(energy_match.group(1))

log_dir = sys.argv[1]
outtbl_path = sys.argv[2]

rows = []
for infile in glob(f'{log_dir}/*.log'):
    mol_id, ext = splitext(basename(infile))
    log_text = open(infile).read()
    energy = log2energy(log_text)
    rows.append({'mol_id': mol_id, 'energy': energy})

pd.DataFrame(rows).to_csv(outtbl_path, index=False)
