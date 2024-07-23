FROM python:3.9

WORKDIR /usr/src/app

# Copy Python scripts and bash script to the container
COPY . /usr/src/app

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# # Make the bash script executable (if needed)
# RUN chmod +x snapshot.sh

# # Expose port if your Python application listens on a specific port
# EXPOSE 5000

# # Define the command to run your Python application
# CMD ["python", "app.py"]

