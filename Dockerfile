FROM python:3.9

WORKDIR /app

# Copy Python scripts and bash script to the container
COPY . /app

# Install Python dependencies
RUN pip install -r requirements.txt

# Make the bash script executable (if needed)
RUN chmod +x snapshot.sh

# Expose port if your Python application listens on a specific port
EXPOSE 8000

# Define the command to run your Python application
CMD ["python", "src/mutantautomate/app.py"]

