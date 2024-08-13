FROM continuumio/miniconda3

# Set the working directory in the container
WORKDIR /usr/src/app

# Copy Python scripts and bash script to the container
COPY . /usr/src/app

# Install pdbfixer using conda
RUN conda install --yes -c conda-forge pdbfixer

# Install build tools
RUN apt-get update && apt-get install -y \
  build-essential \
  gcc \
  g++

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# # Make the bash script executable (if needed)
# RUN chmod +x snapshot.sh

# # Expose port if your Python application listens on a specific port
# EXPOSE 5000

# # Define the command to run your Python application
# CMD ["python", "app.py"]

