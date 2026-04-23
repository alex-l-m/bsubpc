run_csd_python_api search.py
python cif_tblite_forces_to_extxyz.py csd_crystals crystal_forces
python atoms2table.py crystal_forces/*.xyz
Rscript force_histograms.R
