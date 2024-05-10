bg_color white

load /Users/asameerpradhan/MutantAutoMate/src/mutantautomate/5OJ6.pdb

# Color the whole protein green
color green

# Select the residue and color it red
select resn ASP
color red, resn ASP

# Cartoon for residue and sticks for protein
cartoon oval, resn ASP
show ribbon

# Turn the structure for better visibility (rotate around [1,1,1] axis)
turn [1,1,1], 90

set label_size, -3
set label_color, black
set label_font_id, 10
pseudoatom foo
set label_position,(50, 10, 10)
label foo, '5OJ6.pdb'
set ray_opaque_background, 1
ray
png 5OJ6-ASP.png, width=800, height=600, dpi=300
quit
