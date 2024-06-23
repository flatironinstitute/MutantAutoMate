import subprocess
from flask import Flask, request, render_template
from celery import Celery

# Initialize Flask app
app = Flask(__name__)

# Configure Celery
app.config['CELERY_BROKER_URL'] = 'redis://localhost:6379/0'
app.config['CELERY_RESULT_BACKEND'] = 'redis://localhost:6379/0'
celery = Celery(app.name, broker=app.config['CELERY_BROKER_URL'])
celery.conf.update(app.config)

# Define a Celery task
@celery.task
def execute_command(gene_name, residue1, position, residue2, top_isoforms):
    command = f'python3 src/mutantautomate/combined_single.py --gene-name {gene_name} --residue1 {residue1} --position {position} --residue2 {residue2} --top-isoforms {top_isoforms}'
    try:
        process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        stdout, stderr = process.communicate()

        if process.returncode != 0:
            return {'success': False, 'output': stderr}

        return {'success': True, 'output': stdout}

    except subprocess.CalledProcessError as e:
        return {'success': False, 'output': str(e)}

    except Exception as e:
        return {'success': False, 'output': str(e)}

# Flask route using Celery task
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        gene_name = request.form['gene-name']
        residue1 = request.form['residue1']
        position = request.form['position']
        residue2 = request.form['residue2']
        top_isoforms = request.form['top-isoforms']

        # Enqueue Celery task instead of executing command directly
        task = execute_command.delay(gene_name, residue1, position, residue2, top_isoforms)
        return render_template('processing.html', task_id=task.id)

    return render_template('index.html')

# Error handling
@app.errorhandler(404)
def page_not_found(e):
    return render_template('error.html', error='Page not found'), 404

# Run the application
if __name__ == '__main__':
    app.run(debug=True)
