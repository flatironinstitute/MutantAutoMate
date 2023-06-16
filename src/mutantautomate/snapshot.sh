#!/bin/bash

# Get file paths from user input
file1="$1"

# Extract the name of the PDB file
pdb_file=$(basename "$file1")

# Create the script
echo "bg_color white

load $file1

# Select the residue and color it differently
select resi 100
color red, resi 100

# Zoom into the selected residue
zoom resi 100

set label_size, -3
set label_color, black
set label_font_id, 10
pseudoatom foo
set label_position,(50, 10, 10)
label foo, 'TRY'
set ray_opaque_background, 1
ray
png 3.png, width=800, height=600, dpi=300
quit" > try1.pml

# Run the script with PyMOL
pymol -c try1.pml
