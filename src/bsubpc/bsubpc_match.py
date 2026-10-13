"""Shared BsubPc substructure definition and match-index labels."""

TEMPLATE_SMARTS_RAW = (
    "[#5]12-[#7]3:[#6]4:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6]:3-"
    "[#7]=[#6]3:[#6]5:[#6](:[#6]:[#6]:[#6]:[#6]:5):[#6](:[#7]:3-1)="
    "[#7]-[#6]1=[#7]~2-[#6](-[#6]2:[#6]-1:[#6]:[#6]:[#6]:[#6]:2)=[#7]-4"
)

TEMPLATE_SMARTS = TEMPLATE_SMARTS_RAW.replace(":", "~").replace("=", "~").replace("-", "~")

TEMPLATE_ATOM_LABELS = {
    0: "boron",
    1: "pyrrole_nitrogen_1",
    19: "pyrrole_nitrogen_2",
    22: "pyrrole_nitrogen_3",
    10: "imine_nitrogen_1",
    20: "imine_nitrogen_2",
    30: "imine_nitrogen_3",
    6: "outer_terminal_carbon_1",
    7: "outer_terminal_carbon_2",
    15: "outer_terminal_carbon_3",
    16: "outer_terminal_carbon_4",
    27: "outer_terminal_carbon_5",
    28: "outer_terminal_carbon_6",
}

TEMPLATE_ATOM_INDEX = {label: index for index, label in TEMPLATE_ATOM_LABELS.items()}
