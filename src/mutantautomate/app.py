# import os
# import subprocess
# import sys
# from flask import Flask, request

# app = Flask(__name__)

# def create_pml_script(file_path, residue_name, position):
#     # Your existing create_pml_script function code here
#     pass

# def run_pymol_script(script_content):
#     # Your existing run_pymol_script function code here
#     pass

# def generate_snapshot(pdb_file_path, residue_name, position):
#     script_content = create_pml_script(pdb_file_path, residue_name, position)
#     run_pymol_script(script_content)

# @app.route('/')
# def index():
#     gene_name = request.args.get('gene-name')
#     residue1 = request.args.get('residue1')
#     position = request.args.get('position')
#     residue2 = request.args.get('residue2')
#     top_isoforms = request.args.get('top-isoforms')

#     # Assuming combined_single.py is in the same directory as app.py
#     command = f'python3 combined_single.py --gene-name {gene_name} --residue1 {residue1} --position {position} --residue2 {residue2} --top-isoforms {top_isoforms}'

#     # Run the command and get the output
#     result = subprocess.run(command, shell=True, capture_output=True, text=True)

#     # Assuming the PDB file path is passed as an argument to app.py
#     pdb_file_path = sys.argv[1]

#     # Generate snapshot using the PDB file path
#     generate_snapshot(pdb_file_path, residue1, position)

#     return result.stdout

# if __name__ == '__main__':
#     if len(sys.argv) != 2:
#         print("Usage: python app.py /path/to/pdb/file.pdb")
#         sys.exit(1)
    
#     app.run(debug=True)

from flask import Flask

app = Flask(__name__)

@app.route('/')
def index():
    return 'Hello, World!'

if __name__ == '__main__':
    app.run(debug=True)

