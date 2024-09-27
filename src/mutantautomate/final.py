import argparse
import tempfile
from io import BytesIO, StringIO
from pathlib import Path
import json
import re
import urllib
import urllib.parse
import urllib.request
import sys
import subprocess
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import mdtraj as md
from fpdf import FPDF
from Bio import SeqIO, BiopythonDeprecationWarning
from Bio.Align import PairwiseAligner
from pypdb import *
from bioservices import *

import requests
from requests.adapters import HTTPAdapter, Retry
from bs4 import BeautifulSoup
import ast

from PIL import Image as PILImage, ImageFilter
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib import colors, utils
from reportlab.pdfgen import canvas
from reportlab.platypus import (
    SimpleDocTemplate, Image, Spacer, Table, Paragraph, PageTemplate, HRFlowable, ListFlowable, ListItem, Frame
)

# Suppress warnings
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=BiopythonDeprecationWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

# Define regular expression pattern for retrieving next link
re_next_link = re.compile(r'<(.+)>; rel="next"')

# Define retry settings for HTTP requests
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])

# Create a session with retry settings
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


# 1. Function to retrieve the next link from HTTP response headers
def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


# Function for initial search for gene name in Uniprot & retrieve all isoforms for a given gene name
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
if len(sys.argv) != 11:
    print("Usage: python file-name.py --gene-name gene-name --residue1 residue1 --position position --residue2 residue2 --top-isoforms True/False")
    sys.exit(1)

# Create argument parser
parser = argparse.ArgumentParser(description='Description of your program')

# Add arguments
parser.add_argument('--gene-name', type=str, help='Gene name')
parser.add_argument('--residue1', type=str, help='Residue 1')
parser.add_argument('--position', type=int, help='Position')
parser.add_argument('--residue2', type=str, help='Residue 2')
parser.add_argument('--top-isoforms', type=str, help='Show top isoforms: True/False')

# Parse the arguments
args = parser.parse_args()

# Function to retrieve the top 5 isoforms of the gene UniProt ID
def get_top_isoforms(gene_name):
    isoforms = get_all_isoforms(gene_name)
    top_5_isoforms = isoforms[:5]
    top_isoform_list = [isoform for isoform, _ in top_5_isoforms]
    return top_isoform_list

# Call the function to retrieve the top isoforms conditionally iff --top isoforms flag is used
if args.top_isoforms.lower() == 'true':
    top_isoforms = get_top_isoforms(args.gene_name)
    print("Top 5 isoforms of gene UniProt ID:")
    for isoform in top_isoforms:
        print(isoform)

# Call the print_top_isoforms function conditionally
if args.top_isoforms:
    print(get_top_isoforms(args.gene_name))

# Retrieve command-line arguments
gene_name = args.gene_name
residue1 = args.residue1
position = args.position
residue2 = args.residue2
top_isoforms = args.top_isoforms


# Search for matching isoforms
matching_isoforms = search_residue(residue1, position, gene_name)

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
    isoform_sequence = next((isoform[1] for isoform in matching_isoforms if isoform[0] == isoform_id), None)
    isoform_url = f"https://www.uniprot.org/uniprot/{isoform_id}"
else:
    print("Didn't find any matching isoforms for given selection.")
    sys.exit(0)  # Exit the script if no matching isoforms are found


# # Function to calculate similarity between two sequences using PairwiseAligner
# def calculate_similarity(sequence1, sequence2):
#     aligner = PairwiseAligner()
#     aligner.mode = 'global'
#     aligner.match_score = 1
#     aligner.mismatch_score = -1
#     alignments = aligner.align(sequence1, sequence2)
#     best_alignment = alignments[0]
#     alignment_score = best_alignment.score
#     alignment_length = len(best_alignment)
#     similarity = alignment_score / alignment_length * 100
#     return similarity


# Function to score isoforms based on similarity to the gene sequence
# def score_isoforms_by_similarity(gene_name, isoforms):
#     scored_isoforms = []
#     for isoform in isoforms:
#         u = UniProt()
#         entry = u.retrieve(isoform, "fasta")
#         sequence = ''.join(entry.strip().split('\n')[1:])
#         similarity = calculate_similarity(gene_name, sequence)
#         scored_isoforms.append((isoform, similarity))
#     scored_isoforms.sort(key=lambda x: x[1], reverse=True)
#     return scored_isoforms[:3]  # Return the top 3 scored isoforms

# # Retrieve the sequence for the selected isoforms and return the top 3 isoforms
# top_scored_isoforms = score_isoforms_by_similarity(gene_name, matching_isoforms)

# Create UniProt object
u = UniProt()

# Retrieve the sequence for the selected isoform
sequence = u.retrieve(matching_isoforms[0:6], "fasta")
fasta_string = sequence
only_element = matching_isoforms[0]
uniprot_id = only_element[0:6]

# Define the URL for UniProt data
url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"

# Make an HTTP GET request to the URL
response = requests.get(url)

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

else:
    print(f"Error: Failed to retrieve data. Status code: {response.status_code}")

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

# Get the directory of the current file: what to do for a web server???!!!
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

pdbpath = download_pdb(pdb_ids[0], current_dir)

      
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

# Load the trajectory data and topology using MDTraj
traj = md.load(pdbpath)
topology = traj.topology

# Compute the RMSD of each frame in the trajectory to a reference structure
reference_frame = traj[0]  # Assuming the first frame is your reference
rmsd = md.rmsd(traj, reference_frame)


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



# Calculate DSSP assignments
dssp = md.compute_dssp(traj)

# Example: Create a list of structured frames indices
structured_frames_indices = []
for frame_idx, frame_dssp in enumerate(dssp):
    if any(residue_dssp in {'E', 'H'} for residue_dssp in frame_dssp):
        structured_frames_indices.append(frame_idx)

# Example: Check if the target residue is present in the structured frames
target_residue_index = 0  # Set your target residue index here
if any(frame_idx == target_residue_index for frame_idx in structured_frames_indices):
    structured_or_not = f"Residue {position} is in a structured part of the protein"
else:
    structured_or_not = f"Residue {position} is not in a structured part of the protein"

# Generate the final output message
output_message = (
    "This is a mutation from " + str(score) + ".\n" 
    
)

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

    score = calculate_grantham_score(grantham_dict, amino_acid1, amino_acid2)

    threshold = 100  # Define the threshold value for high Grantham score

    if score is not None:
        print(f"The Grantham score between {amino_acids.get(residue1)} and {amino_acids.get(residue2)} is {score}")
        grantham_score= score
        grantham_output = f"The Grantham score between {amino_acids.get(residue1)} and {amino_acids.get(residue2)} is {score}"

        if score > threshold:
            print("This is a high Grantham score, indicating a potentially significant evolutionary distance.")
            grantham_output_extra = "This is a high Grantham score, indicating a potentially significant evolutionary distance."
        else:
            grantham_output_extra = "Not a potentially high score."


# Define the path to your downloaded PDB file
pdb_file = pdbpath  # Assuming `pdbpath` is defined earlier in your script

# Define output HTML file path
output_html = os.path.join(current_dir, "structure.html")

# Run minimal_molview to generate the HTML file
try:
    subprocess.run(["minimal_molview", pdb_file, "-o", output_html], check=True)
    print(f"HTML file generated successfully: {output_html}")
except subprocess.CalledProcessError as e:
    print(f"Error generating HTML file: {e}")
