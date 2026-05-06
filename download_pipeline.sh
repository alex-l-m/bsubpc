set -euo pipefail

# CSD search → mol2 + cif per refcode + match-index tables.
run_csd_python_api search.py
# CSD diagram metadata (used by process_search_results.py to gate cif outputs).
run_csd_python_api extract_csd_diagram_and_molecule_tables.py search_results.csv csd_tables
# RDKit canonicalization (mol output) and ASE crystal-cell read (xyz output).
uv run python process_search_results.py
