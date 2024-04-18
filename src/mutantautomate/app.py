from flask import Flask, request
import subprocess

app = Flask(__name__)

@app.route('/')
def index():
    # Extract arguments from the query string
    gene_name = request.args.get('gene-name')
    residue1 = request.args.get('residue1')
    position = request.args.get('position')
    residue2 = request.args.get('residue2')
    top_isoforms = request.args.get('top-isoforms')

    # Construct the command
    command = f'python3 combined_single.py --gene-name {gene_name} --residue1 {residue1} --position {position} --residue2 {residue2} --top-isoforms {top_isoforms}'

    # Execute the command
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    return result.stdout

if __name__ == '__main__':
    app.run(debug=True)

