FROM continuumio/miniconda3

# Install build tools
RUN apt-get update && apt-get install -y \
  build-essential \
  gcc \
  g++

# Set the working directory in the container
WORKDIR /usr/src/app

# Install Python dependencies
COPY requirements.txt /usr/src/app
RUN pip install --no-cache-dir -r requirements.txt

# Copy Python scripts and bash script to the container
COPY . /usr/src/app

# # Make the bash script executable (if needed)
# RUN chmod +x snapshot.sh

# Expose port if your Python application listens on a specific port
EXPOSE 5000

# Define the command to run your Python application
CMD ["flask", "--app", "src/mutantautomate/app2", "run", "--host=0.0.0.0", "--debug"]

