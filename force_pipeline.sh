uv run python cif_tblite_forces_to_extxyz.py crystal_cells crystal_forces
uv run python check_extxyz_sanity.py
uv run python atoms2table.py crystal_forces/*.xyz
Rscript force_histograms.R
