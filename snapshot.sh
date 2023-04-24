#!/bin/bash


# Get file paths from user input
file1="$1"
file2="$2"

res='131-141'
name='2Z6G'
name2='New_2x6g'
echo ${file1:32:4}
echo $name

#create the script
echo "load $file1
set cartoon_fancy_helices = 1
set cartoon_highlight_color = grey70
bg_color WHITE
set antialias = 1
set ortho = 1
set sphere_mode, 5

colour red, 2Z6G
show cartoon, 2Z6G

select res $res
remove res $res

ray 1000,1500
load $file2
colour red, New_2x6g
align ${file1:32:4}, New_2x6g
set label_size, -3
set label_font_id, 10
pseudoatom foo
set label_position,(3,2,1)
label foo, 'removed part'
png New_2x6g.png" > try1.pml

# GO!
pymol -c try1.pml
