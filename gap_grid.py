"""Render one PNG per BsubPc in csd_matches.csv, sorted by experimental
HOMO-LUMO gap, each split at its boron-axial bond. The companion shell
script gap_grid.sh composites the tiles into a single grid via ImageMagick
montage.

Why this exists: the sTDA model on these crystal geometries mostly just
separates the one outlier with an unusually low gap (XOGYAV, Cl-Cl12BsubPc,
1.94 eV). Putting the structures side by side, sorted by gap, makes the
chemical reason for that outlier visible.

Sources:
  csd_matches.csv     — chemical mol_id ↔ CSD refcode (entry).
  entry_gap.csv       — refcode → experimental gap, output of
                        gap_database.R. Its `mol_id` column is actually
                        the refcode (the R script does
                        `rename(mol_id = entry)`).
  crystal_geom/*.mol  — RDKit-canonicalized structures from the search.
                        Some refcodes in csd_matches.csv aren't here
                        (search didn't find them, or the diagram-atom
                        check in process_search_results.py rejected
                        them); we just skip those rather than fall back,
                        since the outlier is present and the rest of
                        the rows are clustered at 2.14 eV anyway.

Tiles are written as `<rank>_<entry>.png` so a glob in gap_grid.sh
sorts alphabetically (== gap-ascending) and montage lays them out in
the right order without needing a manifest file.
"""
from pathlib import Path

import pandas as pd
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdmolfiles import MolFromMolFile

from diagrams import draw_split_bsubpc
from match import label_core_atoms

TILE_W, TILE_H = 520, 340
TILE_DIR = Path("gap_grid_tiles")

matches = pd.read_csv("csd_matches.csv")
gaps = pd.read_csv("entry_gap.csv").rename(columns={"mol_id": "entry"})
df = matches.merge(gaps, on="entry").sort_values("gap")

TILE_DIR.mkdir(exist_ok=True)
for stale in TILE_DIR.glob("*.png"):
    stale.unlink()

rank = 0
for r in df.itertuples(index=False):
    path = Path(f"crystal_geom/{r.entry}.mol")
    if not path.exists():
        print(f"skip {r.mol_id} ({r.entry}): no mol file")
        continue
    # MolToMolFile drops atom properties, so bsubpc_idx/bsubpc_label have
    # to be recomputed via the SMARTS template after re-reading.
    mol = MolFromMolFile(path, removeHs=False)
    label_core_atoms(mol)

    drawer = rdMolDraw2D.MolDraw2DCairo(TILE_W, TILE_H)
    drawer.drawOptions().legendFontSize = 18
    draw_split_bsubpc(drawer, mol, legend=f"{r.mol_id}  |  {r.gap:.2f} eV  |  {r.entry}")

    tile_path = TILE_DIR / f"{rank:02d}_{r.entry}.png"
    tile_path.write_bytes(drawer.GetDrawingText())
    rank += 1
    print(f"ok   {r.mol_id} ({r.entry}, gap={r.gap:.2f}) -> {tile_path}")

print(f"\nWrote {rank} tiles to {TILE_DIR}/. Run ./gap_grid.sh to composite.")
