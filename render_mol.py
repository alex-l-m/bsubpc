#!/usr/bin/env python3
"""Render a .mol file in PyMOL from the negative Y direction.

Recommended usage (per PyMOL docs):
    pymol -cq render_mol_neg_y.py -- input.mol output.png

This script:
  * loads the input MOL file
  * points the camera from the negative Y side toward the origin
    (so +Z is up and +X is to the viewer's right)
  * fits all atoms in view as tightly as reasonably possible
  * ray-traces to a transparent-background PNG
"""

from __future__ import annotations

import sys
from pathlib import Path

from pymol import cmd


# Easy knobs to tweak later
WIDTH = 1600
HEIGHT = 1200
DPI = 300
ZOOM_BUFFER_ANGSTROM = 0
STICK_RADIUS = 0.18

input_path = Path(sys.argv[1]).expanduser().resolve()
output_path = Path(sys.argv[2]).expanduser().resolve()

if not input_path.exists():
    raise FileNotFoundError(f"Input file does not exist: {input_path}")

output_path.parent.mkdir(parents=True, exist_ok=True)

# Start from a clean object list.
cmd.delete("all")

# Match the working viewport to the final image aspect ratio before zooming.
cmd.viewport(WIDTH, HEIGHT)

# Load the molecule. PyMOL supports MOL input and usually infers the format
# from the file extension.
cmd.load(str(input_path), "mol")

if cmd.count_atoms("mol") == 0:
    raise RuntimeError(f"PyMOL loaded no atoms from: {input_path}")

# Give small molecules a reasonable default look.
cmd.hide("everything", "all")
cmd.show("sticks", "all")
cmd.set("stick_radius", STICK_RADIUS)

# Transparent background for both viewport and ray-traced PNG.
cmd.bg_color("white")
cmd.set("opaque_background", 0)
cmd.set("ray_opaque_background", 0)

# Color carbons
cmd.color("grey50", "elem C")

# Orthoscopic view makes complete=1 fitting more predictable.
cmd.set("orthoscopic", 1)

# Fixed framing knob:
#   smaller value -> zoom in
#   larger value  -> zoom out
FIXED_HALF_SIZE = 6
FIXED_SLAB = 200.0

# Build a temporary, origin-centered reference object.
# We will zoom to THIS object, not to the molecule.
cmd.pseudoatom("_frame", pos=[ FIXED_HALF_SIZE, 0.0, 0.0])
cmd.pseudoatom("_frame", pos=[-FIXED_HALF_SIZE, 0.0, 0.0])
cmd.pseudoatom("_frame", pos=[0.0,  FIXED_HALF_SIZE, 0.0])
cmd.pseudoatom("_frame", pos=[0.0, -FIXED_HALF_SIZE, 0.0])
cmd.pseudoatom("_frame", pos=[0.0, 0.0,  FIXED_HALF_SIZE])
cmd.pseudoatom("_frame", pos=[0.0, 0.0, -FIXED_HALF_SIZE])
cmd.hide("everything", "_frame")

# Reset to the canonical XYZ-aligned view.
cmd.reset()

# Center the window and rotation origin on the literal coordinate origin.
cmd.center("_frame", origin=1)

# Camera on the negative Y side, looking toward the origin.
# Keep your empirically-correct upside-down fix.
cmd.turn("x", 90)
cmd.turn("x", 180)

# Fixed zoom, independent of molecule size.
cmd.zoom("_frame", buffer=0.0, complete=1)

# Give yourself a generous clipping slab so clipping is not the limiting factor.
cmd.clip("slab", FIXED_SLAB)

# Remove the temporary helper object before rendering.
cmd.delete("_frame")

# Ray trace and save.
cmd.png(str(output_path), width=WIDTH, height=HEIGHT, dpi=DPI, ray=1, quiet=1)
