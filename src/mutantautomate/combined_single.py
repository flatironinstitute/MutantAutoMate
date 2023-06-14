import requests
from requests.adapters import HTTPAdapter, Retry
import re
from pypdb import *
from bioservices import *
from pypdb import *
import os
import sys
import urllib
import urllib.request
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
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import SimpleDocTemplate, Paragraph, Image, ListFlowable, ListItem
import subprocess
import warnings

warnings.filterwarnings("ignore", category=UserWarning)


re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_all_isoforms():
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

def search_residue(residue, position):
    matching_isoforms = []
    all_isoforms = get_all_isoforms()
    for isoform, sequence in all_isoforms:
        if len(sequence) > position-1 and sequence[position-1] == residue:
            matching_isoforms.append(isoform)
            #print(isoform)
    return matching_isoforms


input_string = input("Enter the gene name and residue information (e.g., SHANK3 D26): ")
gene_name, residue_info = input_string.split()
residue = residue_info[0]
position = int(residue_info[1:])

matching_isoforms = search_residue(residue, position)
if len(matching_isoforms) > 0:
    #print(f"The residue {residue} is present at position {position} in the following isoform(s):")
    a = "The residue {residue} is present at position {position} in the following isoform(s):"
    for isoform in matching_isoforms:
        random_variable = "I dont know why im getting an err"
else:
    #print(f"The residue {residue} is not present at position {position} in any of the isoforms.")
    a = "The residue {residue} is not present at position {position} in any of the isoforms."


#first uniprot match to fasta sequence
u = UniProt()
sequence = u.retrieve(matching_isoforms[0],"fasta")


from io import StringIO
from Bio import SeqIO

fasta_string = sequence


fasta_io = StringIO(fasta_string) 

records = SeqIO.parse(fasta_io, "fasta") 

for rec in records:
    seq_str = str(rec.seq)
    

fasta_io.close() 


q = Query(seq_str[0:54], 
          query_type="sequence", 
          return_type="polymer_entity")
result = q.search()
if result is None or result is None:
    pdbnewcode = matching_isoforms[0]

else:
    highest_score = -1.0
    identifier = ""
    for result in result['result_set']:
        if result['score'] > highest_score:
            highest_score = result['score']
            identifier = result['identifier']
    #print("Identifier with the highest score:", identifier[0:4])

def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        #print(url)
        #print(outfnm)
        return outfnm
    except Exception as err:
        print("ERROR")
        return  outfnm

# Get the directory of the current file
current_dir = os.path.dirname(os.path.abspath(__file__))

def new_method_for_alphafold(pdbcode, datadir):
    new_url = "https://alphafold.ebi.ac.uk/files/AF-" + matching_isoforms[0] + "-F1-model_v4.pdb"
    pdbfn2 = pdbcode + ".pdb"
    outfnm2 = os.path.join(datadir, pdbfn2)
    #print(new_url)
    try: 
        urllib.request.urlretrieve(new_url, outfnm2)
        #print(outfnm2)
        return outfnm2
    except Exception as err:
        print("ERROR")
        return outfnm2

if result is None or result is None:
    pdbpath = new_method_for_alphafold(matching_isoforms[0], current_dir) #relative path
else:
    pdbpath = download_pdb(identifier[0:4], current_dir)

#getPDF.py
residue2 = input("Enter residue two:")

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
    
   
    #print("Charge change from", f"{residue1_charge}", "to", f"{residue2_charge}")


score = get_charge_change_score(residue, residue2)


# Get the directory of the current file
current_dir2 = os.path.dirname(os.path.abspath(__file__))
#print(current_dir2)

# Define the bash script command
bash_script = os.path.join(current_dir2, "snapshot.sh")

# Extract the filename from the bash_script path
filename = os.path.basename(bash_script)
# Append "snapshot.sh" to the filename
snapshot_script = os.path.join(current_dir2, filename, "snapshot.sh")

file1 = pdbpath

# Run the bash script using subprocess with arguments
output = subprocess.run(["bash", bash_script, file1], capture_output=True, text=True)

# Extract the output path from the command's standard output
screenshot = os.path.join(current_dir2, "3.png")
#print(screenshot)

def generate_pdf(image_path, screenshot_path):
    # Create a new PDF document with letter size
    doc = SimpleDocTemplate("output.pdf", pagesize=letter, rightMargin=50)

    # Define the print statements to be written to the PDF
    print_statements = [f"For the gene {gene_name} the residue at position {position} goes from {residue} to {residue2}.",
                        f"The Uniprot ID for the matched isoform is {isoform}",
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
