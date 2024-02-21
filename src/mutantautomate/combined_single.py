# Import required libraries
from reportlab.lib.pagesizes import letter
from jinja2 import Environment, FileSystemLoader
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Image, Spacer, Table, Paragraph, PageTemplate
from PIL import Image as PILImage, ImageFilter
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Frame
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Spacer
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Image, Paragraph, HRFlowable, ListFlowable, ListItem, Spacer, Table
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import letter
from PIL import Image as PILImage
import tempfile
from reportlab.lib import colors
from io import BytesIO
from reportlab.lib.utils import ImageReader
from reportlab.platypus.flowables import HRFlowable
from reportlab.lib.units import inch
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from PIL import Image as PILImage
from reportlab.lib import utils
from reportlab.platypus.frames import Frame
#remove redundant ones from above^^^ clean code
import requests
from bs4 import BeautifulSoup
import urllib.parse
import json
from requests.adapters import HTTPAdapter, Retry
import re
import ast
from io import StringIO
from Bio import SeqIO
from pypdb import *
from reportlab.platypus import HRFlowable
from bioservices import *
import numpy as np
import os
from Bio import BiopythonDeprecationWarning
import sys
import urllib
import urllib.request
from pathlib import Path
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
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Image, ListFlowable, ListItem, Spacer, PageTemplate
import subprocess
import warnings

from Bio.Align import PairwiseAligner

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)

# Define regular expression pattern for retrieving next link
re_next_link = re.compile(r'<(.+)>; rel="next"')

# Define retry settings for HTTP requests
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])

# Create a session with retry settings
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


# Function to retrieve the next link from HTTP response headers
def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


# Function to retrieve all isoforms for a given gene name
def get_all_isoforms(gene_name):
    isoforms = []
    url1 = "https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+"
    url3 = "&includeIsoform=true&format=list&(taxonomy_id:9606)"
    url = url1 + gene_name + url3
    batch_url = url
    while batch_url:
        response = session.get(batch_url)
        response.raise_for_status()
        for isoform in response.text.strip().split("\n"):
            isoform_url = f"https://www.uniprot.org/uniprot/{isoform}.fasta"
            isoform_response = session.get(isoform_url)
            isoform_response.raise_for_status()
            sequence = ''.join(isoform_response.text.strip().split('\n')[1:])
            isoforms.append((isoform, sequence))
        batch_url = get_next_link(response.headers)
    return isoforms


# Function to search for a specific residue at a position in isoforms of a gene
def search_residue(residue, position, gene_name):
    matching_isoforms = []
    all_isoforms = get_all_isoforms(gene_name)
    for isoform, sequence in all_isoforms:
        if len(sequence) > position - 1 and sequence[position - 1] == residue:
            matching_isoforms.append(isoform)
    return matching_isoforms

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 5:
    print("Usage: python file-name.py gene-name residue1 position residue2")
    sys.exit(1)

# Retrieve command-line arguments
gene_name = sys.argv[1]
residue1 = sys.argv[2]
position = int(sys.argv[3])
residue2 = sys.argv[4]

# Display user inputs
print("Gene Name:", gene_name)
print("Residue 1:", residue1)
print("Position:", position)
print("Residue 2:", residue2)

# Search for matching isoforms
matching_isoforms = search_residue(residue1, position, gene_name)
#print(matching_isoforms)

def get_gene_name(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)
    lines = response.text.split("\n")
    for line in lines:
        if line.startswith("GN   Name="):
            gene_name = line.split("GN   Name=")[1].split(";")[0]
            return gene_name
    return None

# Function to calculate the gene name for each isoform and filter isoforms with different gene names
def filter_isoforms_by_gene(matching_isoforms, gene_name):
    filtered_isoforms = []
    for isoform in matching_isoforms:
        gene_name_isoform = get_gene_name(isoform)
        if gene_name_isoform == gene_name:
            filtered_isoforms.append(isoform)
    return filtered_isoforms

# Filter isoforms by gene name
matching_isoforms = filter_isoforms_by_gene(matching_isoforms, gene_name)

if len(matching_isoforms) > 0:
    # Select the correct isoform
    selected_isoform = max(matching_isoforms, key=lambda isoform: len(isoform))
    isoform_id = selected_isoform[0:6]
    #print(isoform_id)
    isoform_sequence = next((isoform[1] for isoform in matching_isoforms if isoform[0] == isoform_id), None)
    isoform_url = f"https://www.uniprot.org/uniprot/{isoform_id}"
    #print(isoform_url)
else:
    print("Didn't find any matching isoforms for given selection.")
    sys.exit(0)  # Exit the script if no matching isoforms are found

# Rest of the code...

# Function to calculate similarity between two sequences using PairwiseAligner
def calculate_similarity(sequence1, sequence2):
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1
    aligner.mismatch_score = -1
    alignments = aligner.align(sequence1, sequence2)
    best_alignment = alignments[0]
    alignment_score = best_alignment.score
    alignment_length = len(best_alignment)
    similarity = alignment_score / alignment_length * 100
    return similarity


# Function to score isoforms based on similarity to the gene sequence
def score_isoforms_by_similarity(gene_name, isoforms):
    scored_isoforms = []
    for isoform in isoforms:
        u = UniProt()
        entry = u.retrieve(isoform, "fasta")
        sequence = ''.join(entry.strip().split('\n')[1:])
        similarity = calculate_similarity(gene_name, sequence)
        scored_isoforms.append((isoform, similarity))
    scored_isoforms.sort(key=lambda x: x[1], reverse=True)
    return scored_isoforms[:3]  # Return the top 3 scored isoforms

# Retrieve the sequence for the selected isoforms and return the top 3 isoforms
top_scored_isoforms = score_isoforms_by_similarity(gene_name, matching_isoforms)
# Create UniProt object
u = UniProt()

# Retrieve the sequence for the selected isoform
sequence = u.retrieve(matching_isoforms[0:6], "fasta")
fasta_string = sequence
only_element = matching_isoforms[0]
uniprot_id = only_element[0:6]
print(uniprot_id)



# Define the URL for UniProt data
url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"

# Make an HTTP GET request to the URL
response = requests.get(url)
# Rest of the code...

# Check if the request was successful (status code 200)
if response.status_code == 200:
    # Get the response content
    response_text = response.text

    # Function to extract PDB IDs
    def extract_pdb_ids(text):
        # Define a regular expression pattern to match "database":"PDB","id":"<>"
        pattern = r'"database":"PDB","id":"([^"]+)"'

        # Use regex to find all matches of the pattern
        matches = re.findall(pattern, text)

        return matches

    # Extract PDB IDs from the response text
    pdb_ids = extract_pdb_ids(response_text)

    # Print the extracted PDB IDs
    print(pdb_ids[0])
        #print(f"PDB ID: {pdb_id}")
else:
    print(f"Error: Failed to retrieve data. Status code: {response.status_code}")

# Rest of the code...


# Function to download a PDB file from the Internet
def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print("ERROR")
        return outfnm


# Get the directory of the current file
current_dir = os.path.dirname(os.path.abspath(__file__))


# Function to download a PDB file using the new method for AlphaFold
def new_method_for_alphafold(pdbcode, datadir):
    new_url = "https://alphafold.ebi.ac.uk/files/AF-" + uniprot_id + "-F1-model_v4.pdb"
    pdbfn2 = pdbcode + ".pdb"
    outfnm2 = os.path.join(datadir, pdbfn2)
    try:
        urllib.request.urlretrieve(new_url, outfnm2)
        return outfnm2
    except Exception as err:
        print("ERROR")
        return outfnm2


# # Download or use the new method to get the PDB file
# if result is None or result is None:
#     pdbpath = new_method_for_alphafold(uniprot_id, current_dir)  # relative path

pdbpath = download_pdb(pdb_ids[0], current_dir)

# Get user input for residue two
# residue2 = input("Enter residue two:")

# Define the dictionary of amino acid names
amino_acids = {
    'A': 'Alanine',
    'C': 'Cysteine',
    'D': 'Aspartic Acid',
    'E': 'Glutamic Acid',
    'F': 'Phenylalanine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'K': 'Lysine',
    'L': 'Leucine',
    'M': 'Methionine',
    'N': 'Asparagine',
    'P': 'Proline',
    'Q': 'Glutamine',
    'R': 'Arginine',
    'S': 'Serine',
    'T': 'Threonine',
    'V': 'Valine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine'
}

# Function to determine the charge change between two residues
def charge_statement(residue1, residue2):
    aa_charge_dict = {
        "A": "non-polar",
        "C": "polar",
        "D": "negative",
        "E": "negative",
        "F": "bulky",  # "non-polar",
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
        "W": "bulky",  # "non-polar",
        "Y": "bulky"  # "polar"
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

    score = (
        f"{amino_acids.get(residue1)} at position {position} to {amino_acids.get(residue2)}"
        + " from a " + f"{residue1_charge}" + " charged amino acid to a " + f"{residue2_charge}" + " amino acid"
    )
    return score


# Get the charge change score
score = charge_statement(residue1, residue2)

import os
import subprocess
import mdtraj as md

# Define the path for the snapshot script
bash_script = os.path.join(current_dir, "snapshot.sh")
file1 = pdbpath

# Convert one-letter code to three-letter code
amino_acid_mapping = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}

# Convert residue1 one-letter code to three-letter code
residue1_three_letter = amino_acid_mapping.get(residue1, residue1)
# Convert the position to string
converted_position = str(position)

# Run the bash script using subprocess with arguments
output = subprocess.run(["bash", bash_script, pdbpath, residue1_three_letter, converted_position], capture_output=True, text=True)

# Extract the output path from the command's standard output
screenshot = os.path.join(current_dir, "snap.png")
# print(screenshot)

# Load the protein trajectory and topology using mdtraj
traj = md.load(pdbpath)
topology = traj.topology

# Compute the RMSD of each frame in the trajectory to a reference structure
reference_frame = traj[0]  # Assuming the first frame is your reference
rmsd = md.rmsd(traj, reference_frame)


# Compute the SASA for each residue in each frame
sasa = md.shrake_rupley(traj)


# Choose thresholds for RMSD and SASA to determine the structured region
rmsd_threshold = 0.3  # Adjust as per your requirements
sasa_threshold = 10.0  # Adjust as per your requirements


# Ensure that rmsd and sasa have the same number of frames
num_frames = min(rmsd.shape[0], sasa.shape[0])
rmsd = rmsd[:num_frames]
sasa = sasa[:num_frames]


# Find the frames where the RMSD and SASA values are below the thresholds
structured_frames = (rmsd < rmsd_threshold) & (sasa[:, :, np.newaxis] > sasa_threshold)


# Find the residues that are present in the structured frames
structured_residues = []
for residue in topology.residues:
    residue_frames = structured_frames[:, residue.index]
    if residue_frames.any():
        structured_residues.append(residue)

# Define the dictionary of amino acid names
amino_acids = {
    'A': 'Alanine',
    'C': 'Cysteine',
    'D': 'Aspartic Acid',
    'E': 'Glutamic Acid',
    'F': 'Phenylalanine',
    'G': 'Glycine',
    'H': 'Histidine',
    'I': 'Isoleucine',
    'K': 'Lysine',
    'L': 'Leucine',
    'M': 'Methionine',
    'N': 'Asparagine',
    'P': 'Proline',
    'Q': 'Glutamine',
    'R': 'Arginine',
    'S': 'Serine',
    'T': 'Threonine',
    'V': 'Valine',
    'W': 'Tryptophan',
    'Y': 'Tyrosine'
}

target_residue_name = topology.residue(position).name

if target_residue_name is None:
    print(f'Invalid amino acid code: {residue1}')
else:
    # Check if the target residue is present in the structured frames
    target_residue_index = topology.residue(position).index
    residue_frames = structured_frames[:, target_residue_index]
    if residue_frames.any():
        structured_or_not = f"Residue {amino_acids.get(residue1)} is in a structured part of the protein"
    else:
        structured_or_not = f"Residue {amino_acids.get(residue2)} is not in a structured part of the protein"


# Generate the final output message
output_message = (
    "This is a mutation from " + str(score) + ".\n" 
    
)
print(output_message)


def calculate_grantham_score(grantham_dict, aa1, aa2):
    key = (aa1, aa2)
    if key in grantham_dict:
        return grantham_dict[key]
    else:
        print(f"Grantham score not available for ({aa1}, {aa2})")
        return None

if __name__ == "__main__":
    # Load the Grantham dictionary from the output file
    output_file = "grantham_output.txt"
    with open(output_file, "r") as f:
        grantham_dict = ast.literal_eval(f.read())

    # Take user input for amino acids
    amino_acid1 = residue1 #input("Enter the first amino acid: ").upper()
    amino_acid2 = residue2 #input("Enter the second amino acid: ").upper()
    print(amino_acid1, amino_acid2)

    score = calculate_grantham_score(grantham_dict, amino_acid1, amino_acid2)

    threshold = 100  # Define the threshold value for high Grantham score

#Add to PDF
    if score is not None:
        print(f"The Grantham score between {amino_acids.get(residue1)} and {amino_acids.get(residue2)} is {score}")
        grantham_score= score
        grantham_output = f"The Grantham score between {amino_acids.get(residue1)} and {amino_acids.get(residue2)} is {score}"

        if score > threshold:
            print("This is a high Grantham score, indicating a potentially significant evolutionary distance.")
            grantham_output_extra = "This is a high Grantham score, indicating a potentially significant evolutionary distance."
        else:
            grantham_output_extra = "Not a potentially high score."


from weasyprint import HTML

def generate_pdf(image_path, normal_image_path, zoomed_image_path, image_a_path, image_b_path, bg_image_path):
    # HTML content with placeholders for variables
    name = gene_name + "_" + residue1 + str(residue2)
    pdf_filename = f"{name}.pdf"
    html_content = f'''
    <!DOCTYPE html>

    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>MutantAutoMate Report</title>
        <link href="report.css" rel="stylesheet">
        <style>
            body {{
                font-family: 'Fira Sans', sans-serif;
                background-color: #f0f0f0;
                color: #333;
                padding: 20px;
            }}

            header {{
                background-color: #fbc847; /* Update header background color */
                color: #333; /* Update header text color */
                padding: 10px;
                text-align: center;
            }}

            header {{
            background-color: #fbc847;
            color: #333;
            padding: 10px;
            text-align: center;
            position: relative; /* Add relative positioning */
        }}

        .logo {{
            position: absolute; /* Position the logo absolutely within the header */
            top: 10px;
            left: 20px;
            width: 100px; /* Set the width of the logo */
        }}

            nav ul {{
                list-style: none;
                margin: 0;
                padding: 0;
                background-color: #333;
                overflow: hidden;
                margin-bottom: 10px; 
            }}

            nav li {{
          
                margin-bottom: 10px;
            }}

            nav li a {{
                display: block;
                color: white;
                text-align: center;
                padding: 14px 16px;
                text-decoration: none;
            }}

            main {{
                padding: 20px;
            }}

            footer {{
                background-color: #fbc847; /* Update footer background color */
                color: #333; /* Update footer text color */
                padding: 10px;
                text-align: center;
                position: relative; 
            }}
            .image-container {{
                position: relative;
                text-align: center;
                display: flex;
                justify-content: center;
                align-items: center;
                margin: 0 20px; 
            }}

            .image-container img {{
                width: 500px;
                height: 500px;
                display: block;
                margin-bottom: 5px;
                border: 2px solid #333; 
            }}
             article {{
                background-color: #fafafa;
                padding: 10px;
                margin-bottom: 40px;
             }}
             .bullet-list li {{
            background: url("{image_path}") no-repeat left center; /* Set image as bullet point */
            background-size: 50px 50px; /* Set the size of the bullet image */
            padding-left: 70px; /* Adjust padding to make space for the bullet image */
            margin-bottom: 30px; /* Add margin between list items */
        }}
           .bullet-list {{
            list-style: none; /* Remove default list styling */
            padding: 0;
        }}
            #visualization-title,
            #interpretation-title {{
                page-break-before: always;
            }}
            .footer-images {{
            position: absolute;
            bottom: 10px;
            left: 50%;
            transform: translateX(-50%);
            display: flex;
            justify-content: space-between;
            width: 100%;
        }}
        .footer-images img {{
            width: 60px; /* Set the width of the small images */
            height: 30px; /* Set the height of the small images */
            margin: 0 10px; /* Add margin between images */
        }}

            h2 {{
                margin-top: 40px; /* Increase margin for better spacing */
            }}

            h3 {{
                margin-top: 20px; /* Increase margin for better spacing */
            }}
        </style>
    </head>
<body>
    
    <header>
        <h1>MutantAutoMate Report</h1>
    </header>
    


    <main>

            <!-- Add your images here -->
            <div style="text-align: center;">
                <div style="display: inline-block; text-align: center;">
                    <img src="{normal_image_path}" alt="Image A" style="width: 300px; height: 300px; border: 1px solid black; margin-bottom: 10px; display: block;">
                    <p><i>Residue highlighted in red</i></p>
                </div>
                <div style="display: inline-block; text-align: center;">
                    <img src="{zoomed_image_path}" alt="Image B" style="width: 100px; height: 100px; border: 1px solid black; display: block;">
                    <p><i>Zoomed in residue views</i></p>
                </div>
                <br>
                
            </div>
        

            
            <ul>
                <ul class="bullet-list">
                <li>MutantAutoMate is a tool that analyzes mutations in genes.</li>
                <li>The PDB ID and the UniProt ID for the selected isoform are <b>{pdb_ids[0]}</b> and <b>{matching_isoforms[0]}</b>.</li>
                <li>{grantham_output_extra}</li>
                <li>{output_message}</li>
                <li>{structured_or_not}.</li>
        </ul>
            </ul>
            
        </section>
    </main>

    <footer>
          <div class="footer-images">
            <img src="{image_a_path}" alt="Image A Footer">
            <img src="{image_b_path}" alt="Image B Footer">
        </div>
        <p>&copy; 2024 MutantAutoMate. All rights reserved.</p>
    </footer>
</body>
</html>
    '''

    # Generate PDF using WeasyPrint
    HTML(string=html_content, base_url=requests.compat.urljoin('file:', os.getcwd())).write_pdf(pdf_filename, presentational_hints=True)

    print(f"Your PDF report for this mutant has been created: {pdf_filename}")

    # Return the generated PDF filename
    return pdf_filename


# Define a dictionary mapping single-letter amino acid codes to three-letter codes
amino_acid_mapping = {
    'A': 'ALA',
    'R': 'ARG',
    'N': 'ASN',
    'D': 'ASP',
    'C': 'CYS',
    'Q': 'GLN',
    'E': 'GLU',
    'G': 'GLY',
    'H': 'HIS',
    'I': 'ILE',
    'L': 'LEU',
    'K': 'LYS',
    'M': 'MET',
    'F': 'PHE',
    'P': 'PRO',
    'S': 'SER',
    'T': 'THR',
    'W': 'TRP',
    'Y': 'TYR',
    'V': 'VAL',
}

# Get the directory of the current file
current_dir = os.path.dirname(os.path.abspath(__file__))

# Set the image filename
image_filename = "logo.png"
image_aA_path = "FI.png"
image_bB_path = "CCB.png"
zoomed_image_path = "zoom.png"  
bg_image_pathh = "bg_image.png"

# Get the residue code based on the single-letter code

residue_code = amino_acid_mapping.get(residue1, residue1)

# Run the bash script using subprocess with arguments
output = subprocess.run(["bash", bash_script, file1, residue1], capture_output=True, text=True)
#print(output)
# Join the directory path with the image filename
image_path = os.path.join(current_dir, image_filename)
image_a_path = os.path.join(current_dir, image_aA_path)
image_b_path = os.path.join(current_dir, image_bB_path)
normal_image_path = os.path.join(current_dir, f"{pdb_ids[0]}-{residue_code}-{position}.png")
zoomed_image_path = os.path.join(current_dir, f"{pdb_ids[0]}-{residue_code}-{position}-zoom.png")
bg_image_path = os.path.join(current_dir, f"{pdb_ids[0]}-bg_image.png")
print(bg_image_path)
# Call the function to generate the PDF with all the information
generate_pdf(image_path, normal_image_path, zoomed_image_path, image_a_path, image_b_path, bg_image_path)
