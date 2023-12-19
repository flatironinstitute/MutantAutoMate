#!/bin/bash

# Get file path and residue name from user input
file1="$1"
residue_name="$2"

# Convert three-letter residue code to one-letter code
case $residue_name in
    'A') residue_code='ALA' ;;
    'R') residue_code='ARG' ;;
    'N') residue_code='ASN' ;;
    'D') residue_code='ASP' ;;
    'C') residue_code='CYS' ;;
    'Q') residue_code='GLN' ;;
    'E') residue_code='GLU' ;;
    'G') residue_code='GLY' ;;
    'H') residue_code='HIS' ;;
    'I') residue_code='ILE' ;;
    'L') residue_code='LEU' ;;
    'K') residue_code='LYS' ;;
    'M') residue_code='MET' ;;
    'F') residue_code='PHE' ;;
    'P') residue_code='PRO' ;;
    'S') residue_code='SER' ;;
    'T') residue_code='THR' ;;
    'W') residue_code='TRP' ;;
    'Y') residue_code='TYR' ;;
    'V') residue_code='VAL' ;;
    *) residue_code=$residue_name ;;  # Use the entered code if not found in the case statement
esac

# Extract the name of the PDB file
pdb_file=$(basename "$file1")
echo "Residue name: $residue_name"

# Create the script for the normal screenshot
echo "bg_color white

load $file1

# Color the whole protein green
color green

# Select the residue and color it red
select resn ${residue_code}
color red, resn ${residue_code}

# Cartoon for residue and sticks for protein
cartoon oval, resn ${residue_code}
show ribbon

# Turn the structure for better visibility (rotate around [1,1,1] axis)
turn [1,1,1], 90

set label_size, -3
set label_color, black
set label_font_id, 10
pseudoatom foo
set label_position,(50, 10, 10)
label foo, '${pdb_file: (-8)}'
set ray_opaque_background, 1
ray
png ${pdb_file%.*}-${residue_code}.png, width=800, height=600, dpi=300
quit" > try1.pml

# Create the script for the zoomed-in screenshot
echo "bg_color white

load $file1

# Color the whole protein green
color green

# Zoom in on the selected residue
zoom resn ${residue_code}

# Show only the selected residue and color it red
show sticks, resn ${residue_code}
color red, resn ${residue_code}

# Label the zoomed-in residue
set label_size, 20  # Adjusted label size
set label_color, black
set label_font_id, 11
pseudoatom foo
set label_position,(-1, 5, 2)  # Adjusted position
label foo, '${residue_code}'

set ray_opaque_background, 1
ray
png ${pdb_file%.*}-${residue_code}.png, width=800, height=600, dpi=300
png ${pdb_file%.*}-${residue_code}-zoom.png, width=800, height=600, dpi=300
quit" > try1_zoomed.pml

# Run the scripts with PyMOL
pymol -c try1.pml
pymol -c try1_zoomed.pml
