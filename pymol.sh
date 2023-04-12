#!/bin/bash



#create the script
echo "fetch 6KYK
hide everything
set cartoon_fancy_helices = 1
set cartoon_highlight_color = grey70
set antialias = 1
set ortho = 1
set sphere_mode, 5
colour blue, 6kyk
colour red, sampleA
show cartoon, 6kyk
select res 314
ray 1000,1500
zoom res 314
png Example_Output.png" > try1.pml

# GO!
pymol try1.pml