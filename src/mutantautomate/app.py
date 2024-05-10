import os
from flask import Flask, request
import subprocess

app = Flask(__name__)

def create_pml_script(file_path, residue_name):
    residue_code = residue_name_to_code(residue_name)
    pdb_file = os.path.basename(file_path)
    script_content = f"""\
bg_color white

load {file_path}

# Color the whole protein green
color green

# Select the residue and color it red
select resn {residue_code}
color red, resn {residue_code}

# Cartoon for residue and sticks for protein
cartoon oval, resn {residue_code}
show ribbon

# Turn the structure for better visibility (rotate around [1,1,1] axis)
turn [1,1,1], 90

set label_size, -3
set label_color, black
set label_font_id, 10
pseudoatom foo
set label_position,(50, 10, 10)
label foo, '{pdb_file[-8:]}'
set ray_opaque_background, 1
ray
png {pdb_file[:-4]}-{residue_code}.png, width=800, height=600, dpi=300
quit
"""
    return script_content

def residue_name_to_code(residue_name):
    residue_codes = {
        'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
        'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
        'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
    }
    return residue_codes.get(residue_name.upper(), residue_name)

def run_pymol_script(script_content):
    with open('try1.pml', 'w') as f:
        f.write(script_content)
    os.system('pymol -c try1.pml')

def generate_snapshot(file_path, residue_name):
    script_content = create_pml_script(file_path, residue_name)
    run_pymol_script(script_content)

@app.route('/')
def index():
    gene_name = request.args.get('gene-name')
    residue1 = request.args.get('residue1')
    position = request.args.get('position')
    residue2 = request.args.get('residue2')
    top_isoforms = request.args.get('top-isoforms')

    command = f'python3 combined_single.py --gene-name {gene_name} --residue1 {residue1} --position {position} --residue2 {residue2} --top-isoforms {top_isoforms}'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    pdb_file_path = 'path/to/your/file.pdb'  # Update with your actual path
    generate_snapshot(pdb_file_path, residue1)
    
    return result.stdout

if __name__ == '__main__':
    app.run(debug=True)
