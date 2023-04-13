#!/bin/usr/env bash

file1=~/"Downloads/7psu.pdb"
file2=~/"Downloads/New_7bxu.pdb"

res='131-141'
name='7psu'
name2 = 'New_7bxu'
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

colour blue, 7psu
show cartoon, 7psu

select res $res
remove res $res

ray 1000,1500
load $file2
colour red, New_7bxu
align ${file1:32:4}, New_7bxu
zoom res 131-141
png 7psu_try.png" > try1.pml

# GO!
pymol try1.pml