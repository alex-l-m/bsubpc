#!/usr/bin/env python3
"""Compute TBLite forces for CIF files and write force-annotated extxyz files."""

import argparse
from pathlib import Path

from ase.io import read, write
from tblite.ase import TBLite


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cif_files = input_dir.glob("*.cif")
    if not cif_files:
        raise SystemExit(f"No CIF files found in {input_dir}")

    for cif_path in cif_files:
        atoms = read(cif_path)
        atoms.set_pbc((True, True, True))
        atoms.calc = TBLite(
            method="GFN2-xTB",
            charge=0,
            multiplicity=1,
            verbosity=0,
        )

        atoms.get_forces()

        out_path = output_dir / f"{cif_path.stem}.xyz"
        write(
            out_path,
            atoms,
            format="extxyz",
            columns=["symbols", "positions", "forces"],
            write_results=True,
        )
        print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
