from Bio import PDB

def mutate_residue(pdb_filename, position, original_residue, new_residue, output_filename):
    # Create a PDB parser object
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_filename)
    
    # Create a PDB IO object to save the modified structure
    io = PDB.PDBIO()
    
    # Function to mutate residues
    def mutate_residue_callback(structure):
        for model in structure:
            for chain in model:
                for res in chain:
                    if res.get_id()[1] == position and res.get_resname() == original_residue:
                        # Create a new residue with the new residue name
                        res.resname = new_residue
                        print(f"Mutating residue at position {position} from {original_residue} to {new_residue}")
    
    # Call the mutation function
    mutate_residue_callback(structure)
    
    # Save the modified structure to a new PDB file
    io.set_structure(structure)
    io.save(output_filename)

# Define the parameters
pdb_filename = '3biw.pdb'         # Input PDB file
position = 140                     # Position of the residue to mutate
original_residue = 'ASP'           # Original residue name (Aspartic Acid)
new_residue = 'TYR'                # New residue name (Tyrosine)
output_filename = 'mutated.pdb'    # Output PDB file

# Run the mutation
mutate_residue(pdb_filename, position, original_residue, new_residue, output_filename)

