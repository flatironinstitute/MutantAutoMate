from fpdf import FPDF
from matplotlib.backends.backend_pdf import PdfPages
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import mdtraj as md
from matplotlib import pyplot as plt
import pandas as pd

residue1 = input("Enter residue one:")
residue2 = input("Enter residue two:")
number = int(input("Enter residue position:"))

#charge change
#1 charge
def get_charge_change_score(residue1, residue2):
    aa_charge_dict = {
        "A": "non-polar",
        "C": "polar",
        "D": "negative",
        "E": "negative",
        "F": "bulky",#"non-polar",
        "G": "non-polar",
        "H": "positive",
        "I": "non-polar",
        "K": "positive",
        "L": "non-polar",
        "M": "non-polar",
        "N": "polar",
        "P": "non-polar",
        "Q": "polar",
        "R": "positive",
        "S": "polar",
        "T": "polar",
        "V": "non-polar",
        "W": "bulky",#"non-polar",
        "Y": "bulky"#"polar"
    }
    
    aa_charge_categories = {
        "positive-to-negative": ["K", "R", "H"],
        "positive-to-hydrophobic": ["K", "R", "H"],
        "negative-to-positive": ["D", "E"],
        "negative-to-hydrophobic": ["D", "E"],
        "hydrophobic-to-positive": ["A", "F", "G", "I", "L", "M", "P", "V", "W", "Y"],
        "hydrophobic-to-negative": ["A", "F", "G", "I", "L", "M", "P", "V", "W", "Y"],
        "hydrophobic-to-polar": ["A", "F", "G", "I", "L", "M", "P", "V", "W", "Y"],
        "polar-to-hydrophobic": ["C", "N", "Q", "S", "T", "Y"],
        "polar-to-positive": ["C", "N", "Q", "S", "T", "Y"],
        "polar-to-negative": ["C", "N", "Q", "S", "T", "Y"],
        "non-polar-to-polar": ["A", "F", "G", "I", "L", "M", "P", "V", "W", "Y"]

    }
    
    residue1_charge = aa_charge_dict[residue1]
    residue2_charge = aa_charge_dict[residue2]
    
    # print(f"residue1: {residue1}")
    # print(f"residue2: {residue2}")
    print("Charge change from", f"{residue1_charge}", "to", f"{residue2_charge}")


score = get_charge_change_score(residue1, residue2)
#print(score)
#DSSP secondary structure

import matplotlib.pyplot as plt
import mdtraj as md

import seaborn as sns
sns.set(style='white')
sns.set_context('talk')

traj_file = '/Users/asameerpradhan/Downloads/6kyk.pdb'
top_file = '/Users/asameerpradhan/Downloads/6kyk.pdb'

traj = md.load(traj_file,top=top_file)

topology = traj.topology

dssp = md.compute_dssp(traj)

if dssp[0][number] == "H":
    print("alpha-helix structure")
else:
    print("not an alpha helix")



# #SASA graphs(x2) snapshot
traj = md.load("/Users/asameerpradhan/Downloads/6kyk.pdb")
sasa = md.shrake_rupley(traj, mode='atom')
print(type(plt.plot(sasa[0])))

# #VDW
# #image construction

# #add 2 PDF

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

# Define the function to generate the PDF
def generate_pdf():
    # Create a new PDF document with letter size
    c = canvas.Canvas("output.pdf", pagesize=letter)

    # Define the print statements to be written to the PDF
    print_statements = [f"the residue changes from {score}" , "This is a sample PDF generated with Python."]

    # Set the starting position for writing the print statements
    x, y = 50, 750

    # Write the print statements to the PDF
    for statement in print_statements:
        c.drawString(x, y, statement)
        y -= 50

    # Save the PDF document
    c.save()

# Call the function to generate the PDF
generate_pdf()

import mdtraj as md

# Load your PDB file
traj1 = md.load("/Users/asameerpradhan/Downloads/6kyk.pdb") #main PDB
traj2 = md.load("/Users/asameerpradhan/Downloads/7sk2.pdb") #mutated PDB from getPDB code
# Compute the SASA for each frame
sasa1 = md.shrake_rupley(traj1, mode='atom')
sasa2 = md.shrake_rupley(traj2, mode='atom')
print(type(sasa[0]))

#generate SASA PDF

data_1 = sasa1[0]
data_2 = sasa2[0]

df_1 = pd.DataFrame(data_1)
df_2 = pd.DataFrame(data_2)

with PdfPages(r'/Users/asameerpradhan/output1.pdf') as export_pdf:
    plt.plot(sasa1[0])
    plt.xlabel('residue', fontsize=14)
    plt.ylabel('sasa', fontsize=14)
    plt.grid(True)
    export_pdf.savefig()
    plt.close()

with PdfPages(r'/Users/asameerpradhan/output2.pdf') as export_pdf:
    plt.plot(sasa2[0])
    plt.xlabel('residue(mutated)', fontsize=14)
    plt.ylabel('sasa(mutated)', fontsize=14)
    plt.grid(True)
    export_pdf.savefig()
    plt.close()

