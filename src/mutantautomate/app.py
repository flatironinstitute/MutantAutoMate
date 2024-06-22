import subprocess
from flask import Flask, request, render_template

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Retrieve form data
        gene_name = request.form['gene-name']
        residue1 = request.form['residue1']
        position = request.form['position']
        residue2 = request.form['residue2']
        top_isoforms = request.form['top-isoforms']
        
        # Construct command to execute
        command = f'python3 combined_single.py --gene-name {gene_name} --residue1 {residue1} --position {position} --residue2 {residue2} --top-isoforms {top_isoforms}'
        
        try:
            # Run the command and capture output
            result = subprocess.run(command, shell=True, capture_output=True, text=True, stderr=subprocess.STDOUT)
            
            # Check for return code to handle errors
            if result.returncode != 0:
                error_message = f"Error executing command: {result.stdout}"
                return render_template('error.html', error=error_message)
            
            # Assuming success, render the result template
            return render_template('result.html', response=result.stdout)
        
        except subprocess.CalledProcessError as e:
            # Handle subprocess error (command execution error)
            error_message = f"Error executing command: {e}"
            return render_template('error.html', error=error_message)
        
        except Exception as e:
            # Handle any other unexpected errors
            error_message = f"Unexpected error: {e}"
            return render_template('error.html', error=error_message)
    
    # Render the HTML form template for GET requests
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
