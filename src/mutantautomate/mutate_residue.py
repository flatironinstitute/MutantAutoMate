from Bio import PDB
from io import StringIO


def mutate_residue(pdb_string, position, original_residue, new_residue):
    # Create a PDB parser object
    parser = PDB.PDBParser(QUIET=True)

    # Create a filehandle from the string
    pdb_handle = StringIO(pdb_string)

    return "hahahaha"
