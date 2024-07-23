import json
import re
import time
from flask import (
    Flask,
    render_template,
    request,
    jsonify,
    Response,
    stream_with_context,
)
import requests
from requests.adapters import HTTPAdapter, Retry


app = Flask(__name__)

# Define retry settings for HTTP requests
# Create a session with retry settings
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session = requests.Session()
session.mount("https://", HTTPAdapter(max_retries=retries))


@app.route("/")
def index():
    return render_template("index2.html")


@app.route("/stream-test")
def stream_test():
    def generate():
        page = 1
        while True:
            results = get_paginated_results(page)
            if not results:
                print("No more results")
                break
            yield f"data: {json.dumps(results)}\n\n"
            page += 1
            time.sleep(1)  # Simulate delay between pages

    return Response(stream_with_context(generate()), content_type="text/event-stream")


def get_paginated_results(page):
    # Simulate a database query that returns paginated results
    # This is just an example; replace with your actual data fetching logic
    data = [
        {"id": 1, "name": "Item 1"},
        {"id": 2, "name": "Item 2"},
        {"id": 3, "name": "Item 3"},
        # Add more items as needed
    ]
    page_size = 1
    start = (page - 1) * page_size
    end = start + page_size
    return data[start:end]


@app.route("/search-uniprot", methods=["GET"])
def search_uniprot_route():
    gene_name = request.args.get("gene_name")
    result = search_uniprot(gene_name)
    return jsonify(result)


def search_uniprot(gene_name):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=reviewed:true+AND+{gene_name}&includeIsoform=true&format=list&(taxonomy_id:9606)"
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
    return {"isoforms": isoforms, "next_link": next_link}


@app.route("/grantham", methods=["GET"])
def grantham():
    import os
    import ast

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

    error = None

    # Get values of "amino_acid_1" and "amino_acid_2" from the query string
    amino_acid_1 = request.args.get("amino_acid_1")
    amino_acid_2 = request.args.get("amino_acid_2")

    # Get the path of the current file's directory
    current_directory = os.path.dirname(__file__)
    # Get the path of "grantham_output.txt" in the current directory
    grantham_output_path = os.path.join(current_directory, "grantham_output.txt")
    # Load the Grantham dictionary from the output file
    with open(grantham_output_path, "r") as f:
        grantham_dict = ast.literal_eval(f.read())

    # Calculate the Grantham score
    grantham_score = None
    grantham_key = (amino_acid_1, amino_acid_2)
    if grantham_key in grantham_dict:
        grantham_score = grantham_dict[grantham_key]
    else:
        error = f"Grantham score not found for {amino_acid_1} and {amino_acid_2}"

    threshold = 100  # Define the threshold value for high Grantham score

    grantham_output = None
    grantham_output_extra = None

    if grantham_score is not None:
        grantham_output = f"The Grantham score between {amino_acids.get(amino_acid_1)} and {amino_acids.get(amino_acid_2)} is {grantham_score}"

        if grantham_score > threshold:
            grantham_output_extra = "This is a high Grantham score, indicating a potentially significant evolutionary distance."
        else:
            grantham_output_extra = "Not a potentially high score."

    # Return the values as JSON
    return jsonify(
        {
            "amino_acid_1": amino_acid_1,
            "amino_acid_2": amino_acid_2,
            "grantham_score": grantham_score,
            "grantham_output": grantham_output,
            "grantham_output_extra": grantham_output_extra,
            "error": error,
        }
    )
