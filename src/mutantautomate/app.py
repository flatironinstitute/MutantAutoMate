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
        command = f'python3 src/mutantautomate/combined_single.py --gene-name {gene_name} --residue1 {residue1} --position {position} --residue2 {residue2} --top-isoforms {top_isoforms}'
        
        try:
            # Run the command and capture output
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate()

            # Check for return code to handle errors
            if process.returncode != 0:
                error_message = f"Error executing command: {stderr}"
                return render_template('error.html', error=error_message)
            
            # Assuming success, render the result template
            return render_template('result.html', response=stdout)
        
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
