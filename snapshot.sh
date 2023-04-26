#!/bin/bash


# Get file paths from user input
file1="$1"
file2="$2"
echo $file2
echo ${file2:32:5}

#create the script
echo "bg_color white  

load $file1
colour blue, ${file1:32:4}

load $file2
color red, ${file2:32:5}

align ${file1:32:4}, ${file2:32:5}

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

# GO!
pymol -c try1.pml
