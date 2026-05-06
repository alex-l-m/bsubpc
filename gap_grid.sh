#!/usr/bin/env bash
# Composite the per-molecule tiles written by gap_grid.py into one PNG.
#
# Tiles are named gap_grid_tiles/<rank>_<entry>.png so the alphabetical
# glob is in gap-ascending order, matching the iteration order of
# gap_grid.py.
#
# Options:
#   -tile 4x          : 4 columns; rows auto-computed.
#   -geometry +5+5    : keep tiles at native size, with 5 px of horizontal
#                       and vertical spacing on each side (i.e. 10 px
#                       between adjacent tiles, 5 px outer border).
#   -background white : fill spacing and any uneven trailing cells.
#   -font DejaVu-Sans : the only font on this Nix-managed VM that
#                       imagemagick auto-detects; without an explicit
#                       -font, montage's text subsystem prints a noisy
#                       (non-fatal) "unable to read font" warning even
#                       though no labels are drawn.

set -euo pipefail

magick montage \
    gap_grid_tiles/*.png \
    -tile 4x \
    -geometry +5+5 \
    -background white \
    -font DejaVu-Sans \
    gap_ordered_grid.png

echo "Wrote gap_ordered_grid.png from $(ls gap_grid_tiles/*.png | wc -l) tiles"
