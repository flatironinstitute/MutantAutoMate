from flask import Flask, render_template

app = Flask(__name__)


@app.route("/")
def index():
    return render_template("index2.html")


@app.route("/grantham", methods=["GET"])
def grantham():
    from flask import request, jsonify
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

    print("grantham dict", grantham_dict)

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
            print(
                "This is a high Grantham score, indicating a potentially significant evolutionary distance."
            )
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

    # threshold = 100  # Define the threshold value for high Grantham score

    # if score is not None:
    #     print(
    #         f"The Grantham score between {amino_acids.get(residue1)} and {amino_acids.get(residue2)} is {score}"
    #     )
    #     grantham_score = score
    #     grantham_output = f"The Grantham score between {amino_acids.get(residue1)} and {amino_acids.get(residue2)} is {score}"

    #     if score > threshold:
    #         print(
    #             "This is a high Grantham score, indicating a potentially significant evolutionary distance."
    #         )
    #         grantham_output_extra = "This is a high Grantham score, indicating a potentially significant evolutionary distance."
    #     else:
    #         grantham_output_extra = "Not a potentially high score."
    # # Return grantham_score, grantham_output, grantham_output_extra as JSON
    # return jsonify(
    #     {
    #         "grantham_score": grantham_score,
    #         "grantham_output": grantham_output,
    #         "grantham_output_extra": grantham_output_extra,
    #     }
    # )
