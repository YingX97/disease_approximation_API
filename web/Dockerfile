# Use the official Python image from the Docker Hub
FROM python:3.11-slim

# Set environment variables to avoid interactive prompts during package installation
ENV PYTHONUNBUFFERED=1

# Set the working directory in the container
WORKDIR /app

# Copy the requirements.txt file into the container at /app
COPY requirements.txt /app/

# Install the dependencies
RUN pip install -r requirements.txt

# Copy the rest of the application code into the container at /app
COPY . /app

# Command to run the Flask application
CMD gunicorn -b :8080 app:app
