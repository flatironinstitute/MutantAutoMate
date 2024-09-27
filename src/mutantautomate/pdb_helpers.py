from Bio import PDB
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO
from io import StringIO
import pdbfixer
from openmm.app import PDBFile
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

def mutate_residue_2(pdb_string, chain_id, position, to_residue_one_letter):
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


def fixer_to_string(fixer):
    output_file_handle = StringIO()
    PDBFile.writeFile(
        fixer.topology, fixer.positions, file=output_file_handle, keepIds=True
    )
    out_string = output_file_handle.getvalue()
    return out_string


def mutate_residue(pdb_string, position, original_residue_, new_residue_, chain_id="A"):
    # Convert the residue names to three-letter codes
    original_residue = amino_acid_mapping[original_residue_]
    new_residue = amino_acid_mapping[new_residue_]

    # Create a filehandle from the string
    pdb_file_handle = StringIO(pdb_string)

    # Create a mutation string based on the position and the new residue
    mutation_string = f"{original_residue}-{position}-{new_residue}"

    fixer = pdbfixer.PDBFixer(pdbfile=pdb_file_handle)
    fixer.applyMutations([mutation_string], chain_id)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()

    out_string = fixer_to_string(fixer)

    return out_string


def trim_pdb(pdb_string, chains_to_keep):
    # Create a filehandle from the string
    pdb_file_handle = StringIO(pdb_string)
    fixer = pdbfixer.PDBFixer(pdbfile=pdb_file_handle)
    chains = fixer.topology.chains()
    chains_list = list(chains)
    print(f"Number of chains: {len(chains_list)}")
    chains_to_remove = []

    # print(list(chains), len(list(chains)))

    for chain in chains_list:
        print(chain.id)
        this_chain_id = chain.id
        if this_chain_id not in chains_to_keep:
            print(f"Removing chain {this_chain_id}")
            chains_to_remove.append(this_chain_id)
        else:
            print(f"Keeping chain {this_chain_id}")

    fixer.removeChains(chainIds=chains_to_remove)

    out_string = fixer_to_string(fixer)

    return out_string


# def mutate_residue_2(pdb_string):
#     output_file_handle = StringIO(pdb_string)
#     parser = PDBParser(QUIET=1)
#     structure = parser.get_structure("my_structure", output_file_handle)
#     # structure = parser.get_structure("my_structure", "5e0m.pdb")
#     # mutate(structure, 'A', 140, 'TYR', mutation_type='first')
#     mutate(structure, 'A', 140, 'TYR')
#     # all_atoms = list(structure.get_atoms())
#     return "haha"