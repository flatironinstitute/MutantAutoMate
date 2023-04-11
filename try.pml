#########################
### Load your protein ###
#########################

load /Users/asameerpradhan/Downloads/6kyk.pdb, 6kyk

##########################
### Set your viewpoint ###
##########################

#################
### Set Style ###
#################

hide everything
set cartoon_fancy_helices = 1
set cartoon_highlight_color = grey70
bg_colour white
set antialias = 1
set ortho = 1
set sphere_mode, 5

############################
### Make your selections ###
############################

select sampleA, 6kyk and resi 314

colour blue, 6kyk
colour red, sampleA
show cartoon, 6kyk


###################
### Save a copy ###
###################

ray 1000,1500
png Lysozyme_Example_Output.png