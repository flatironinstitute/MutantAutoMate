#!/bin/bash

# Get file path and residue number from user input
file1="$1"
residue_number="$2"

# Extract the name of the PDB file
pdb_file=$(basename "$file1")
echo "$residue_number"

# Create the script
echo "bg_color white

load $file1

# Color the whole protein green
color green

# Select the residue and color it red
select resi ${residue_number}
color red, res ${residue_number}

#cartoon for residue and sticks for protein
cartoon oval, res ${residue_number}
show sticks

#rotate for better visibility
rotate [1,1,1], ${residue_number}, chain A

set label_size, -3
set label_color, black
set label_font_id, 10
pseudoatom foo
set label_position,(50, 10, 10)
label foo, '${file1: (-8)}'
set ray_opaque_background, 1
ray
png snap.png, width=800, height=600, dpi=300
quit" > try1.pml

# Run the script with PyMOL
pymol -c try1.pml
