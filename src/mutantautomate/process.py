import os
import re
import ast
import requests
from requests.adapters import HTTPAdapter, Retry

from Bio.Align import PairwiseAligner

# Define retry settings for HTTP requests
# Create a session with retry settings
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def calculate_grantham_score(residue1, residue2):
    # Get the path of the current file's directory
    current_directory = os.path.dirname(__file__)
    # Get the path of "grantham_output.txt" in the current directory
    grantham_output_path = os.path.join(current_directory, "grantham_output.txt")
    # Load the Grantham dictionary from the output file
    with open(grantham_output_path, "r") as f:
        grantham_dict = ast.literal_eval(f.read())
    key = (residue1, residue2)
    if key in grantham_dict:
        return grantham_dict[key]
    else:
        print(f"Grantham score not available for ({residue1}, {residue2})")
        return None


# Define the dictionary of amino acid names
amino_acids = {
    "A": "Alanine",
    "C": "Cysteine",
    "D": "Aspartic Acid",
    "E": "Glutamic Acid",
    "F": "Phenylalanine",
    "G": "Glycine",
    "H": "Histidine",
    "I": "Isoleucine",
    "K": "Lysine",
    "L": "Leucine",
    "M": "Methionine",
    "N": "Asparagine",
    "P": "Proline",
    "Q": "Glutamine",
    "R": "Arginine",
    "S": "Serine",
    "T": "Threonine",
    "V": "Valine",
    "W": "Tryptophan",
    "Y": "Tyrosine",
}


def get_grantham_score_with_statement(residue1, residue2):
    # Calculate the Grantham score
    grantham_score = calculate_grantham_score(residue1, residue2)
    threshold = 100  # Define the threshold value for high Grantham score
    grantham_output = f"Grantham score not available for ({residue1}, {residue2})"
    grantham_output_extra = ""
    if grantham_score is not None:
        grantham_output = f"The Grantham score between {amino_acids.get(residue1)} and {amino_acids.get(residue2)} is {grantham_score}."
        if grantham_score > threshold:
            grantham_output_extra = "This is a high Grantham score, indicating a potentially significant evolutionary distance."
        else:
            grantham_output_extra = "Not a potentially high score."
    # Return the values as JSON
    return {
        "grantham_score": grantham_score,
        "grantham_statement": grantham_output + " " + grantham_output_extra,
    }


# Function to determine the charge change between two residues
def get_charge_statement(residue1, residue2, position):
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
        "Y": "bulky",  # "polar"
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
        "non-polar-to-polar": ["A", "F", "G", "I", "L", "M", "P", "V", "W", "Y"],
    }

    residue1_charge = aa_charge_dict[residue1]
    residue2_charge = aa_charge_dict[residue2]

    charge_statement = (
        f"This is a mutation "
        f"from {amino_acids.get(residue1)} "
        f"at position {position} "
        f"to {amino_acids.get(residue2)}. "
        f"This is a mutation from a {residue1_charge} charged amino acid "
        f"to a {residue2_charge} amino acid."
    )
    return charge_statement


def search_uniprot(url):
    print(f"Searching Uniprot: {url}")
    response = session.get(url)
    response.raise_for_status()
    isoforms = response.text.strip().split("\n")
    next_link = None
    if "Link" in response.headers:
        # Define regular expression pattern for retrieving next link
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        match = re_next_link.match(response.headers["Link"])
        if match:
            next_link = match.group(1)
    return {
        "type": "uniprot_results",
        "isoforms": isoforms,
        "next_link": next_link,
    }


def search_uniprot_generator(gene_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=(reviewed:true)+AND+(taxonomy_id:9606)+AND+{gene_name}&includeIsoform=true&format=list"
    while url:
        result = search_uniprot(url)
        yield result
        url = result["next_link"]


def get_sequence(isoform):
    url = f"https://www.uniprot.org/uniprot/{isoform}.fasta"
    response = session.get(url)
    response.raise_for_status()
    sequence = "".join(response.text.strip().split("\n")[1:])
    return sequence


def get_sequences_generator(isoforms):
    for isoform in isoforms:
        sequence = get_sequence(isoform)
        yield {"type": "isoform_and_sequence", "isoform": isoform, "sequence": sequence}


def search_residue(sequence, residue, position):
    if len(sequence) > position - 1 and sequence[position - 1] == residue:
        return True
    return False


def get_gene_name(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.json"
    response = session.get(url)
    response.raise_for_status()
    data = response.json()
    gene_name = None
    if "genes" in data:
        first_gene = data["genes"][0]
        try:
            gene_name = first_gene["geneName"]["value"]
        except KeyError:
            gene_name = None
    return gene_name


def get_gene_names_generator(isoforms):
    for isoform in isoforms:
        gene_name = get_gene_name(isoform)
        yield {
            "type": "isoform_and_gene_name",
            "isoform": isoform,
            "gene_name": gene_name,
        }


# Function to calculate similarity between two sequences using PairwiseAligner
def calculate_similarity(sequence1, sequence2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    alignments = aligner.align(sequence1, sequence2)
    best_alignment = alignments[0]
    alignment_score = best_alignment.score
    return alignment_score


def process(gene_name, residue1, position, residue2):
    print("Gene Name:", gene_name)
    print("Residue 1:", residue1)
    print("Position:", position)
    print("Residue 2:", residue2)

    # Grantham Score
    grantham_score_with_statement = get_grantham_score_with_statement(
        residue1, residue2
    )
    yield {
        "type": "grantham_score",
        "message": "Calculated Grantham score.",
        **grantham_score_with_statement,
    }

    # Charge Statement
    charge_statement = get_charge_statement(residue1, residue2, position)
    yield {
        "type": "charge_statement",
        "charge_statement": charge_statement,
        "message": "Got charge statement.",
    }

    # Collect all isoforms
    all_isoforms = {}
    for result in search_uniprot_generator(gene_name):
        isoforms = result["isoforms"]
        for isoform in isoforms:
            all_isoforms[isoform] = {
                "sequence": None,
                "gene_name": None,
            }
        yield {
            "type": "isoforms_progress",
            "message": f"Fetched {len(all_isoforms)} isoforms.",
        }
    yield {
        "type": "all_isoforms",
        "message": "Fetched all isoforms.",
        "all_isoforms": list(all_isoforms.keys()),
    }

    # Collect all sequences
    sequences_found = 0
    for result in get_sequences_generator(all_isoforms):
        isoform = result["isoform"]
        sequence = result["sequence"]
        all_isoforms[isoform]["sequence"] = sequence
        sequences_found += 1
        yield {
            "type": "isoform_sequence",
            "message": f"Fetched sequence for {isoform}.",
            "isoform": isoform,
        }
        yield {
            "type": "isoform_sequence_progress",
            "message": f"Fetched {sequences_found} / {len(all_isoforms)} sequences.",
        }

    # Find isoforms with residue1 at position
    matching_isoforms = {}
    for isoform, data in all_isoforms.items():
        sequence = data["sequence"]
        if search_residue(sequence, residue1, position):
            matching_isoforms[isoform] = data
            yield {
                "type": "message",
                "message": f"Found {isoform} with {residue1} at position {position}.",
            }
    yield {
        "type": "matching_isoforms",
        "message": f"Found {len(matching_isoforms)} matching isoforms.",
        "matching_isoforms": list(matching_isoforms.keys()),
    }

    # Get gene names for matching isoforms
    for result in get_gene_names_generator(matching_isoforms.keys()):
        isoform = result["isoform"]
        this_gene_name = result["gene_name"]
        all_isoforms[isoform]["gene_name"] = this_gene_name
        yield {
            "type": "gene_name",
            "message": f"Fetched gene name for {isoform}: {this_gene_name}.",
            "gene_name": this_gene_name,
        }

    # Filter isforms that match the input gene name
    filtered_isoforms = {}
    for isoform in matching_isoforms:
        this_gene_name = matching_isoforms[isoform]["gene_name"]
        yield {
            "type": "message",
            "message": f"Found gene names for {isoform}: {this_gene_name} and {gene_name}.",
        }
        if this_gene_name == gene_name:
            filtered_isoforms[isoform] = matching_isoforms[isoform]
    yield {
        "type": "filtered_isoforms",
        "message": f"Filtered isoforms: {filtered_isoforms.keys()}.",
        "filtered_isoforms": list(filtered_isoforms.keys()),
    }

    # Pairwise similary scores
    pairwise_scores = {}
    for isoform1 in filtered_isoforms:
        for isoform2 in filtered_isoforms:
            sequence1 = filtered_isoforms[isoform1]["sequence"]
            sequence2 = filtered_isoforms[isoform2]["sequence"]
            alignment_score = calculate_similarity(sequence1, sequence2)
            pairwise_scores[(isoform1, isoform2)] = alignment_score
            yield {
                "type": "message",
                "message": f"Computed similarity for {isoform1} and {isoform2}: {alignment_score}.",
            }
    # Convert the dictionary with tuple keys to a list of dictionaries that is JSON serializable
    pairwise_scores_list = []
    for (isoform1, isoform2), score in pairwise_scores.items():
        pairwise_scores_list.append(
            {"isoform1": isoform1, "isoform2": isoform2, "score": score}
        )
    yield {
        "type": "pairwise_scores",
        "message": "Computed pairwise scores.",
        "pairwise_scores": pairwise_scores_list,
    }

    # PDB IDs
    pdb_ids = {}
    for isoform in filtered_isoforms:
        pdb_ids[isoform] = []
        json_url = f"https://www.uniprot.org/uniprot/{isoform}.json"
        response = session.get(json_url)
        response.raise_for_status()
        data = response.json()
        references = data["uniProtKBCrossReferences"]
        for reference in references:
            if reference["database"] == "PDB":
                pdb_id = reference["id"]
                chains_specifier = None
                resolution = None
                if "properties" in reference:
                    for prop in reference["properties"]:
                        if prop["key"] == "Chains":
                            chains_specifier = prop["value"]
                        if prop["key"] == "Resolution":
                            resolution = prop["value"]
                pdb_ids[isoform].append([pdb_id, chains_specifier, resolution])
    yield {"type": "pdb_ids", "message": "Fetched PDB IDs.", "pdb_ids": pdb_ids}
