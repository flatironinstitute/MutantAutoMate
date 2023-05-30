import matplotlib.pyplot as plt
import mdtraj as md
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import PyPDF2


#original ref: '/Users/asameerpradhan/Downloads/6kyk.pdb'
#WT ref: '/Users/asameerpradhan/Downloads/7sk2.pdb'

sns.set(style='white')
sns.set_context('talk')

residue1 = input("Enter residue one:")
residue2 = input("Enter residue two:")
number = int(input("Enter residue position:"))
original_file = input("Enter WT file path: ")
mutated_file = input("Enter mutated file path: ")

traj_file = original_file
top_file = original_file

traj = md.load(traj_file,top=top_file)

topology = traj.topology

dssp = md.compute_dssp(traj)


# #SASA graphs(x2) snapshot
traj = md.load(original_file)
sasa = md.shrake_rupley(traj, mode='atom')
print(type(plt.plot(sasa[0])))
import mdtraj as md

# Load your PDB file
traj1 = md.load(original_file) #main PDB
traj2 = md.load(mutated_file) #mutated PDB from getPDB code
# Compute the SASA for each frame
sasa1 = md.shrake_rupley(traj1, mode='atom')
sasa2 = md.shrake_rupley(traj2, mode='atom')
print(type(sasa[0]))

#generate SASA PDF

data_1 = sasa1[0]
data_2 = sasa2[0]

df_1 = pd.DataFrame(data_1)
df_2 = pd.DataFrame(data_2)

#SASA difference check:
sasa_diff = sasa2[0][number] - sasa1[0][72]
print(sasa_diff)


with PdfPages(r'/Users/asameerpradhan/akshada.pdf') as export_pdf:
    plt.figure(figsize=(3, 3))
    plt.plot(sasa1[0])
    plt.xlabel('residue', fontsize=14)
    plt.ylabel('sasa', fontsize=14)
    plt.grid(True)
    export_pdf.savefig()
    plt.close()
    
    plt.plot(sasa2[0])
    plt.xlabel('residue(mutated)', fontsize=14)
    plt.ylabel('sasa(mutated)', fontsize=14)
    plt.grid(True)
    export_pdf.savefig()
    plt.close()




from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

# Define the function to generate the PDF
def generate_pdf(print_statements):
    # Create a new PDF document with letter size
    c = canvas.Canvas("output.pdf", pagesize=letter)

    # Set the starting position for writing the print statements
    x, y = 50, 750

    # Set the font size for the print statements
    c.setFont("Helvetica", 14)

    # Write the print statements to the PDF
    for statement in print_statements:
        # Split the statement into words
        words = statement.split()

        # Initialize an empty string to store the current line
        current_line = ""

        # Iterate through the words in the statement
        for word in words:
            # Calculate the width of the current line with the new word
            current_line_width = c.stringWidth(current_line + " " + word)

            # If the current line is too long, write it to the PDF and start a new line
            if current_line_width > 500:
                c.drawString(x, y, current_line)
                y -= 25  # Move to the next line
                current_line = ""

            # Add the word to the current line
            current_line += " " + word

        # Write the remaining words in the current line to the PDF
        c.drawString(x, y, current_line)
        y -= 25  # Move to the next line

    # Save the PDF document
    c.save()
if sasa_diff>=0.3:
    # Define the list of print statements
    print_statements = ["The change in SASA is considerable: {sasa_diff}",
                        "Reasoning: The solvent accessible surface area (SASA) score is a measure of the exposed surface area of a molecule that is accessible to solvent molecules.",
                        "It can be used to assess the structural changes that occur in proteins due to mutations and how these changes may impact protein function, stability, and pathogenicity.",
                        "Generally, larger SASA score changes are more likely to have a significant impact on protein function and pathogenicity."]
else: 
    print_statements = ["The change in SASA isn't considerable enough to affect mutant pathogenicity.",
                        "The solvent accessible surface area (SASA) score is a measure of the exposed surface area of a molecule that is accessible to solvent molecules.",
                        "It can be used to assess the structural changes that occur in proteins due to mutations and how these changes may impact protein function, stability, and pathogenicity.",
                        "Generally, larger SASA score changes are more likely to have a significant impact on protein function and pathogenicity."]
# Call the function to generate the PDF
generate_pdf(print_statements)





