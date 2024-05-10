bg_color white

load /Users/asameerpradhan/MutantAutoMate/src/mutantautomate/5OJ6.pdb

# Color the whole protein green
color green

# Zoom in on the selected residue
zoom resn ASP

# Show only the selected residue and color it red
show sticks, resn ASP
color red, resn ASP

# Label the zoomed-in residue
set label_size, 20  # Adjusted label size
set label_color, black
set label_font_id, 11
pseudoatom foo
set label_position,(-1, 5, 2)  # Adjusted position
label foo, 'ASP'

set ray_opaque_background, 1
ray
png 5OJ6-ASP.png, width=800, height=600, dpi=300
png 5OJ6-ASP-zoom.png, width=800, height=600, dpi=300
quit
