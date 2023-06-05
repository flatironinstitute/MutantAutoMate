from fpdf import FPDF
from matplotlib.backends.backend_pdf import PdfPages
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
import mdtraj as md
from matplotlib import pyplot as plt
import pandas as pd
import os
import mdtraj as md
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib import utils
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Image, ListFlowable, ListItem
import subprocess

#user input gene name
gene_name = input("Enter the gene you want to search for (e.g., SHANK3): ")

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


# Get the directory of the current file
current_dir2 = os.path.dirname(os.path.abspath(__file__))
print(current_dir2)

# Define the bash script command
bash_script = os.path.join(current_dir2, "snapshot.sh")

# Extract the filename from the bash_script path
filename = os.path.basename(bash_script)
# Append "snapshot.sh" to the filename
snapshot_script = os.path.join(current_dir2, filename, "snapshot.sh")

file1 = "/Users/asameerpradhan/Downloads/nl_xray_B_D2.pdb"

# Run the bash script using subprocess with arguments
output = subprocess.run(["bash", bash_script, file1], capture_output=True, text=True)

# Extract the output path from the command's standard output
screenshot = os.path.join(current_dir2, "3.png")
print(screenshot)

def generate_pdf(image_path, screenshot_path):
    # Create a new PDF document with letter size
    doc = SimpleDocTemplate("output.pdf", pagesize=letter, rightMargin=50)

    # Define the print statements to be written to the PDF
    print_statements = [f"For the gene {gene_name} the residue at position {number} goes from {residue1} to {residue2}.",
                        f"Parameters that may contribute to the pathogenicity of the mutant are: charge change, presence on alpha-helix strand, and change in solvent accessible surface area.",
                        f"{score}.",
                        f"It's important to note that the specific effect of a point mutation on the pathogenicity of a mutation will depend on many factors, including the nature of the amino acid change, the location of the mutation within the alpha helix, the role of the affected residue in protein function, and the overall protein context.",
                        f"Detailed experimental or computational analysis is typically required to accurately assess the impact of a point mutation on pathogenicity in a specific protein."]

    # Create a list to hold the flowables (Paragraphs, Images, and Bullet List)
    flowables = []

    # Load and add the logo image to the flowables
    logo = utils.ImageReader(image_path)
    logo_width, logo_height = logo.getSize()
    logo_scale = 0.1  # Adjust the scale as needed
    logo_width *= logo_scale
    logo_height *= logo_scale
    img = Image(image_path, width=logo_width, height=logo_height)
    img.hAlign = 'RIGHT'  # Align the image to the right side
    img.top = doc.pagesize[1] - logo_height  # Position the image at the top-right corner
    flowables.append(img)

    # Add the title at the top of the PDF
    styles = getSampleStyleSheet()
    title_text = "MutantAutoMate"
    title = Paragraph(title_text, styles["Title"])
    flowables.append(title)

    # Create a Bullet List flowable
    bullet_list = ListFlowable(
        [
            ListItem(Paragraph(statement, styles["Bullet"])) for statement in print_statements
        ],
        bulletType="bullet",
        leftIndent=40,
        rightIndent=10,
        start='bulletchar',
        bulletColor='black',
        bulletFontName='Helvetica',
        bulletFontSize=10,
        bulletOffsetY=-2,
        bulletDedent='auto',
        spaceAfter=12
    )
    flowables.append(bullet_list)

    # Load and add the screenshot image to the flowables
    screenshot = utils.ImageReader(screenshot_path)
    screenshot_width, screenshot_height = screenshot.getSize()
    screenshot_scale = 0.5  # Adjust the scale as needed
    screenshot_width *= screenshot_scale
    screenshot_height *= screenshot_scale
    screenshot_img = Image(screenshot_path, width=screenshot_width, height=screenshot_height)
    screenshot_img.hAlign = 'CENTER'  # Align the image to the center
    flowables.append(screenshot_img)

    # Build the PDF document with the flowables
    doc.build(flowables)

# Get the directory of the current file
current_dir = os.path.dirname(os.path.abspath(__file__))

# Set the image filename
image_filename = "logo.png"

# Join the directory path with the image filename
image_path = os.path.join(current_dir, image_filename)

# Call the function to generate the PDF
generate_pdf(image_path, screenshot)
