from Bio import PDB
from io import StringIO
import pdbfixer
from openmm.app import PDBFile

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


def mutate_residue(pdb_string, position, original_residue_, new_residue_):
    # Convert the residue names to three-letter codes
    original_residue = amino_acid_mapping[original_residue_]
    new_residue = amino_acid_mapping[new_residue_]

    # Create a filehandle from the string
    pdb_file_handle = StringIO(pdb_string)

    # Create a mutation string based on the position and the new residue
    mutation_string = f"{original_residue}-{position}-{new_residue}"
    chain_id = "A"

    fixer = pdbfixer.PDBFixer(pdbfile=pdb_file_handle)
    fixer.applyMutations([mutation_string], chain_id)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    # fixer.addMissingHydrogens()

    # out_string = pdb_file_handle.getvalue()
    output_file_handle = StringIO()
    PDBFile.writeFile(
        fixer.topology, fixer.positions, file=output_file_handle, keepIds=True
    )

    out_string = output_file_handle.getvalue()

    return out_string


# def mutate_residue(pdb_string, position, original_residue_, new_residue_):
#     # Create a PDB parser object
#     parser = PDB.PDBParser(QUIET=True)

#     # Create a filehandle from the string
#     pdb_handle = StringIO(pdb_string)

#     structure = parser.get_structure("structure", pdb_handle)

#     # Convert the residue names to three-letter codes
#     original_residue = amino_acid_mapping[original_residue_]
#     new_residue = amino_acid_mapping[new_residue_]

#     for model in structure:
#         for chain in model:
#             for res in chain:
#                 if (
#                     res.get_id()[1] == position
#                     and res.get_resname() == original_residue
#                 ):
#                     # Create a new residue with the new residue name
#                     res.resname = new_residue
#                     print(
#                         f"Mutating residue at position {position} from {original_residue} to {new_residue}"
#                     )

#     output_file_handle = StringIO()

#     # Create a PDB IO object to save the modified structure
#     out_structure = PDB.PDBIO()
#     out_structure.set_structure(structure)
#     out_structure.save(output_file_handle)

#     out_string = output_file_handle.getvalue()

#     return out_string
