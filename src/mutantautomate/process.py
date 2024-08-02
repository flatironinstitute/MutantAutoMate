import re
import requests
from requests.adapters import HTTPAdapter, Retry

from Bio.Align import PairwiseAligner

# Define retry settings for HTTP requests
# Create a session with retry settings
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


def construct_uniprot_search_url(gene_name):
    # base_url = "https://rest.uniprot.org/uniprotkb/search"
    # params = {
    #     "query": f"reviewed:true AND {gene_name}",
    #     "includeIsoform": "true",
    #     "format": "list",
    #     "(taxonomy_id:9606)": "",
    # }
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {
        "query": f"reviewed:true AND taxonomy_id:9606 AND {gene_name}",
        "includeIsoform": "true",
        "format": "list",
    }
    full_url = requests.Request("GET", base_url, params=params).prepare().url
    return full_url


def search_uniprot(url):
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
        "isoforms": isoforms,
        "next_link": next_link,
    }


def search_uniprot_generator(gene_name):
    url = construct_uniprot_search_url(gene_name)
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
        yield {"isoform": isoform, "sequence": sequence}


# # Function to search for a specific residue at a position in isoforms of a gene
# def search_residue(all_sequences, residue, position, gene_name):
#     matching_isoforms = []
#     for isoform, sequence in all_sequences:
#         if len(sequence) > position - 1 and sequence[position - 1] == residue:
#             matching_isoforms.append(isoform)
#     return matching_isoforms


# Function to search for a specific residue at a position in isoforms of a gene
# def search_residue(all_isoforms, residue, position, gene_name):
#     matching_isoforms = {}
#     for isoform, data in all_isoforms.items():
#         sequence = data["sequence"]
#         if len(sequence) > position - 1 and sequence[position - 1] == residue:
#             matching_isoforms[isoform] = data
#     return matching_isoforms


def search_residue(sequence, residue, position):
    if len(sequence) > position - 1 and sequence[position - 1] == residue:
        return True
    return False


def get_gene_name(uniprot_id):
    url = f"https://www.uniprot.org/uniprot/{uniprot_id}.txt"
    response = requests.get(url)
    lines = response.text.split("\n")
    for line in lines:
        if line.startswith("GN   Name="):
            gene_name = line.split("GN   Name=")[1].split(";")[0]
            return gene_name
    return None


def get_gene_names_generator(isoforms):
    for isoform in isoforms:
        gene_name = get_gene_name(isoform)
        yield {"isoform": isoform, "gene_name": gene_name}


# Function to calculate similarity between two sequences using PairwiseAligner
def calculate_similarity(sequence1, sequence2):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 1
    aligner.mismatch_score = -1
    alignments = aligner.align(sequence1, sequence2)
    print(f"number of alignments: {len(alignments)}")
    best_alignment = alignments[0]
    # print("best alignment:")
    # print(best_alignment)
    alignment_score = best_alignment.score
    print(f"alignment_score: {alignment_score}")
    alignment_length = len(best_alignment)
    print(f"alignment_length: {alignment_length}")
    print(f"sequence1 length: {len(sequence1)}")
    print(f"sequence2 length: {len(sequence2)}")
    similarity = alignment_score / alignment_length * 100
    return similarity


def process(gene_name, residue1, position, residue2, top_isoforms):
    """Main entry point for the script."""
    print("Call your main application code here")
    print("Gene Name:", gene_name)
    print("Residue 1:", residue1)
    print("Position:", position)
    print("Residue 2:", residue2)
    print("Top Isoforms:", top_isoforms)

    # Collect all isoforms
    all_isoforms = {}
    for result in search_uniprot_generator(gene_name):
        isoforms = result["isoforms"]
        for isoform in isoforms:
            all_isoforms[isoform] = {
                "sequence": None,
                "gene_name": None,
            }
        yield f"got {len(all_isoforms)} isoforms"

    # Collect all sequences
    sequences_found = 0
    for result in get_sequences_generator(all_isoforms):
        isoform = result["isoform"]
        sequence = result["sequence"]
        all_isoforms[isoform]["sequence"] = sequence
        sequences_found += 1
        yield f"got sequence for {isoform}"
        yield f"got {sequences_found} / {len(all_isoforms)} sequences"

    # Find isoforms with residue1 at position
    matching_isoforms = {}
    for isoform, data in all_isoforms.items():
        sequence = data["sequence"]
        if search_residue(sequence, residue1, position):
            matching_isoforms[isoform] = data
            yield f"found {isoform} with {residue1} at position {position}"
    yield f"found {len(matching_isoforms)} matching isoforms: {matching_isoforms.keys()}"

    # Get gene names for matching isoforms
    for result in get_gene_names_generator(matching_isoforms.keys()):
        isoform = result["isoform"]
        gene_name = result["gene_name"]
        all_isoforms[isoform]["gene_name"] = gene_name
        yield f"got gene name for {isoform}: {gene_name}"

    # Filter isforms that match the input gene name
    filtered_isoforms = {}
    for isoform in matching_isoforms:
        print(f"isoform: {isoform}")
        this_gene_name = matching_isoforms[isoform]["gene_name"]
        print(f"this_gene_name: {this_gene_name}")
        if this_gene_name == gene_name:
            print(f"it matches {isoform} {this_gene_name}")
            filtered_isoforms[isoform] = matching_isoforms[isoform]
        else:
            print(f"it doesn't match {isoform} {this_gene_name}")
    yield f"filtered isoforms: {filtered_isoforms.keys()}"

    # Pairwise similary scores
    pairwise_scores = {}
    for isoform1 in filtered_isoforms:
        print(f"isoform1: {isoform1}")
        for isoform2 in filtered_isoforms:
            print(f"isoform2: {isoform2}")
            sequence1 = filtered_isoforms[isoform1]["sequence"]
            sequence2 = filtered_isoforms[isoform2]["sequence"]
            similarity = calculate_similarity(sequence1, sequence2)
            pairwise_scores[(isoform1, isoform2)] = similarity
            yield f"similarity for {isoform1} and {isoform2}: {similarity}"

    print(pairwise_scores)

    # TODO: Is this right?
    # Get the sequence of the main gene?
    # input_gene_sequence = get_sequence(gene_name)
    # print(f"input gene sequence: {input_gene_sequence}")

    # Score isoforms by similarity
    # for isoform, data in filtered_isoforms:
    #     sequence = data["sequence"]
    #     similarity = calculate_similarity(sequence, sequence)
    #     yield f"similarity for {isoform}: {similarity}"

    # # Score isoforms by similarity
    # scored_isoforms = []
    # for isoform in filtered_isoforms:
    #     # Get the sequence from all_sequences
    #     sequence = next(
    #         (seq for isoform_id, seq in all_sequences if isoform_id == isoform), None
    #     )

    # # Test similartity
    # similarty = calculate_similarity(all_sequences[0][1], all_sequences[1][1])
    # yield f"similarity: {similarty}"
    # print("done")


def main():
    results = process("NLGN1", "D", 140, "Y", True)
    for result in results:
        print(result)


main()


# # Collect all sequences
# all_sequences = []
# for result in get_sequences_generator(all_isoforms):
#     isoform = result["isoform"]
#     sequence = result["sequence"]
#     yield f"got sequence for {isoform}"
#     all_isoforms[isoform] = sequence
#     all_sequences.append((isoform, sequence))
#     yield f"got {len(all_sequences)} / {len(all_isoforms)} sequences"

# # Find isoforms with residue1 at position
# matching_isoforms = search_residue(all_sequences, residue1, position, gene_name)
# yield f"found {len(matching_isoforms)} matching isoforms: {matching_isoforms}"

# # Get gene names for isoforms
# isoforms_with_gene_names = []
# for result in get_gene_names_generator(matching_isoforms):
#     print(f"got gene name: {result}")
#     isoforms_with_gene_names.append((result["isoform"], result["gene_name"]))
