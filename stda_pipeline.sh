#!/usr/bin/env bash
# sTDA pipeline: take the per-molecule mol2 files extracted from CSD by
# search.py and produce xtb4stda + stda excitation-energy logs and a single
# summary CSV. Runs after `download_pipeline.sh`.

set -euo pipefail

MOL2_DIR=csd_molecules
XYZ_DIR=crystal_geom_xyz
LOG_DIR=crystal_geom_stda_logs
WORK_DIR=crystal_geom_stda_work
SUMMARY_CSV=crystal_stda.csv.gz

mkdir -p "$XYZ_DIR" "$LOG_DIR" "$WORK_DIR"

# Convert mol2 -> xyz. Restricted to mol_ids that passed the formula gate in
# process_search_results.py (mol_status == 'ok'): the matched component's
# element-count signature must equal one of the CSD chemical-diagram view
# formulas. Same standard the cif side applies via csd_diagram_checks.csv.
uv run python mol2_to_xyz.py \
    "$MOL2_DIR" "$XYZ_DIR" \
    --summary processing_summary.csv

# xtb4stda writes wfn.xtb / charges / energy / gradient into cwd, and stda
# reads wfn.xtb from cwd, so each molecule runs in its own subdirectory.
# Parallelize across molecules; pin each xtb4stda/stda invocation to one
# native thread so total parallelism stays close to the host core count.
WORKERS=${WORKERS:-$(nproc)}

run_one() {
    local INFILE=$1
    local MOLID
    MOLID=$(basename "$INFILE" .xyz)
    local SUBDIR=$WORK_DIR/$MOLID
    rm -rf "$SUBDIR"
    mkdir -p "$SUBDIR"
    cp "$INFILE" "$SUBDIR/$MOLID.xyz"
    if (
        cd "$SUBDIR" && \
        OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 xtb4stda "$MOLID.xyz" > xtb4stda.out 2>&1 && \
        OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 stda -xtb > stda.out 2>&1
    ); then
        cp "$SUBDIR/stda.out" "$LOG_DIR/$MOLID.log"
        printf 'ok   %s\n' "$MOLID"
    else
        printf 'FAIL %s\n' "$MOLID"
    fi
}
export -f run_one
export XYZ_DIR LOG_DIR WORK_DIR

find "$XYZ_DIR" -maxdepth 1 -name '*.xyz' -print0 | \
    xargs -0 -P "$WORKERS" -I '{}' bash -c 'run_one "$@" || true' _ '{}'

n_ok=$(find "$LOG_DIR" -maxdepth 1 -name '*.log' 2>/dev/null | wc -l)
n_total=$(find "$XYZ_DIR" -maxdepth 1 -name '*.xyz' 2>/dev/null | wc -l)
echo "stda: $n_ok / $n_total ok"

uv run python summarize_stda.py "$LOG_DIR" "$SUMMARY_CSV"

# Visualization step. gap_grid.py reads crystal_geom/*.mol files and silently
# skips refcodes whose mol files aren't there, so it tolerates the sTDA set
# being a superset of the canonicalized set.
uv run python gap_database.py
uv run python gap_grid.py
if ls gap_grid_tiles/*.png >/dev/null 2>&1; then
    ./gap_grid.sh
fi
