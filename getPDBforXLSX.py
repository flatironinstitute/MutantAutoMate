import requests
import re
from requests.adapters import HTTPAdapter, Retry
from pathlib import Path
from bioservices import *
from pypdb import *
import os
import urllib.request
from io import StringIO
from Bio import SeqIO
import pandas as pd
from io import StringIO
from Bio import SeqIO
from Bio import ExPASy
from Bio import SwissProt

re_next_link = re.compile(r'<(.+)>; rel="next"')
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))

def get_next_link(headers):
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)

def get_all_isoforms(gene_name): #works!
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

# def search_residue(residue, position, all_isoforms):
#   matching_isoforms = []
#   for i in range(len((all_isoforms))):
#     print(len((all_isoforms)))
#     #iterate through entire list
#     isoform, sequence = all_isoforms[i] #one more loop to check in all sequences
#     #print(isoform, sequence)
#     if len(sequence) > position-1 and sequence[position-1] == residue:
#         matching_isoforms.append(isoform)  
#     return matching_isoforms

def search_residue(residue, position, all_isoforms):
    matching_isoforms = []
    for isoform, sequence in all_isoforms:
        #print(isoform, sequence)
        if len(sequence) > position - 1 and sequence[position - 1] == residue:
            matching_isoforms.append(isoform)
    return matching_isoforms



def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        return  outfnm

def new_method_for_alphafold(pdbcode, datadir):
    new_url = "https://alphafold.ebi.ac.uk/files/AF-" + pdbcode + "-F1-model_v4.pdb"
    pdbfn2 = pdbcode + ".pdb"
    outfnm2 = os.path.join(datadir, pdbfn2)
    try:
        urllib.request.urlretrieve(new_url, outfnm2)
       
        return outfnm2
    except Exception as err:
        return outfnm2

def retrieve_fasta(matching_isoforms):
    sequences = []
    for isoform in matching_isoforms:
        with ExPASy.get_sprot_raw(isoform) as handle:
            record = SwissProt.read(handle)
            sequences.append((isoform, record.sequence))
    return sequences

# Read Excel file
df = pd.read_excel('/Users/asameerpradhan/Downloads/sample.xlsx')

# Extract gene name, residue name, and residue position as separate lists
gene_names = df['gene_name'].tolist()
residue_names = df['residue_name'].tolist()
residue_positions = df['residue_position'].astype(int).tolist()

# Iterate over gene_names list and call the listed methods for each gene name
for (gene_name, residue_name, residue_position, i) in zip(gene_names, residue_names, residue_positions, range(len(residue_names))):
    all_isoforms = get_all_isoforms(gene_name)
    matching_isoforms = search_residue(residue_name, residue_position, all_isoforms)
    #print(matching_isoforms[0])
    #first uniprot match to fasta sequence
    u = UniProt()
    sequence = u.retrieve(matching_isoforms[0],"fasta")
    #print("RETRIEVE ONLY", sequence)
    fasta_string = sequence #select only the sequence part
    fasta_io = StringIO(fasta_string)  #convert to instance of a string object
    records = SeqIO.parse(fasta_io, "fasta")  #iterator object
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
        print("Identifier with the highest score:", identifier[0:4])

     
    if result is None or result is None:
        pdbpath = new_method_for_alphafold(matching_isoforms[0], "/Users/asameerpradhan/Desktop")
    else:
        pdbpath = download_pdb(identifier[0:4], "/Users/asameerpradhan/Desktop")
