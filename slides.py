#!/usr/bin/env python3

import sys
from pathlib import Path


def mol_ids_from_stdin():
    for line in sys.stdin:
        mol_id = line.strip()
        if mol_id:
            yield mol_id


def slide_block(mol_id: str) -> str:
    return "\n".join([
        f"## {mol_id}",
        "",
        "::: {layout-ncol=2}",
        f"![](<diagrams/{mol_id}.svg>)",
        "",
        f"![](<3d_diagrams/{mol_id}.png>)",
        ":::",
    ])


def main() -> None:
    lines = [
        "---",
        "format: revealjs",
        "---",
        "",
    ]

    first_slide = True

    for mol_id in mol_ids_from_stdin():
        diagram_path = Path("diagrams") / f"{mol_id}.svg"
        diagram_3d_path = Path("3d_diagrams") / f"{mol_id}.png"

        if not (diagram_path.is_file() and diagram_3d_path.is_file()):
            continue

        if not first_slide:
            lines.extend([
                "---",
                "",
            ])

        lines.append(slide_block(mol_id))
        lines.append("")
        first_slide = False

    sys.stdout.write("\n".join(lines))


if __name__ == "__main__":
    main()
