from io import StringIO
from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from pymut import mutate

# Convert one-letter code to three-letter code
amino_acid_mapping = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

def mutate_residue(pdb_string, chain_id, position, to_residue_one_letter):
    input_handle = StringIO(pdb_string)
    parser = PDBParser(QUIET=1)
    structure = parser.get_structure("my_structure", input_handle)
    # Mutate the structure
    to_residue_three_letter = amino_acid_mapping[to_residue_one_letter]
    mutate(structure, chain_id, position, to_residue_three_letter)
    # Create a PDBIO object
    io = PDBIO()
    io.set_structure(structure)
    # Save the structure to the StringIO object
    output_handle = StringIO()
    io.save(output_handle)
    # Retrieve the PDB string from the StringIO object
    output_pdb_string = output_handle.getvalue()
    # Return the PDB string
    return output_pdb_string


def trim_pdb(pdb_string, chains_to_keep):

    output_file_handle = StringIO(pdb_string)
    parser = PDBParser(QUIET=1)
    structure = parser.get_structure("my_structure", output_file_handle)

    # Iterate through models and chains
    for model in structure:
        chains_to_remove = []
        for chain in model:
            # Check if the chain ID is in the list of chains to keep
            if chain.id not in chains_to_keep:
                model.detach_child(chain.id)

    # Create an output PDB string
    output_handle = StringIO()
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_handle)
    out_string = output_handle.getvalue()
    
    return out_string
