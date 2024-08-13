from Bio import PDB
from io import StringIO


def mutate_residue(pdb_string, position, original_residue, new_residue):
    # Create a PDB parser object
    parser = PDB.PDBParser(QUIET=True)

    # Create a filehandle from the string
    pdb_handle = StringIO(pdb_string)

    structure = parser.get_structure("structure", pdb_handle)

    for model in structure:
        for chain in model:
            for res in chain:
                if (
                    res.get_id()[1] == position
                    and res.get_resname() == original_residue
                ):
                    # Create a new residue with the new residue name
                    res.resname = new_residue
                    print(
                        f"Mutating residue at position {position} from {original_residue} to {new_residue}"
                    )

    output_file_handle = StringIO()

    # Create a PDB IO object to save the modified structure
    out_structure = PDB.PDBIO()
    out_structure.set_structure(structure)
    out_structure.save(output_file_handle)

    out_string = output_file_handle.getvalue()

    return out_string
