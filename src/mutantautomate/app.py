from flask import Flask, request, render_template
import subprocess

app = Flask(__name__)

# Function to run the combined_single.py script with provided arguments
def run_combined_single(gene_name, residue1, position, residue2, top_isoforms):
    try:
        # Construct command to execute
        command = f'python3 combined_single.py --gene-name {gene_name} --residue1 {residue1} --position {position} --residue2 {residue2} --top-isoforms {top_isoforms}'
        
        # Run the command and capture output
        result = subprocess.run(command, shell=True, capture_output=True, text=True, stderr=subprocess.STDOUT)
        
        # Check for return code to handle errors
        if result.returncode != 0:
            return None, f"Error executing command: {result.stdout}"
        
        # Assuming success, return the result
        return result.stdout, None
    
    except subprocess.CalledProcessError as e:
        return None, f"Error executing command: {e}"
    
    except Exception as e:
        return None, f"Unexpected error: {e}"

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        try:
            # Retrieve form data
            gene_name = request.form['gene-name']
            residue1 = request.form['residue1']
            position = request.form['position']
            residue2 = request.form['residue2']
            top_isoforms = request.form['top-isoforms']
            
            # Call function to run combined_single.py
            output, error = run_combined_single(gene_name, residue1, position, residue2, top_isoforms)
            
            if error:
                return render_template('error.html', error=error)
            
            # Assuming success, render the result template with output
            return render_template('result.html', response=output)
        
        except Exception as e:
            # Handle any unexpected errors
            error_message = f"Unexpected error: {e}"
            return render_template('error.html', error=error_message)
    
    # Render the HTML form template for GET requests
    return render_template('index.html')

if __name__ == '__main__':
    app.run(debug=True)
