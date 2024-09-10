import json
import os
from dotenv import load_dotenv
from google.cloud import storage
from flask import Response, request
from flask_restful import Resource, abort

# helper functions that does the data processing
from models.highest_measurement import (
    compute_gene_measurement,
    get_highest_measurement,
)

from models.files import process_h5_file

from api.v1.exceptions import model_exceptions

# load environment variables from a .env file into the operating system's environment variables.
load_dotenv()


class HighestMeasurement(Resource):
    """
    API Resource class to handle requests for the highest expression measurement
    of a given feature (gene) across multiple all diseases and their datasets . It returns the top n
    disease and cell types combination with the highest expression of the specified gene.
    """

    @model_exceptions
    def post(self):
        args = request.args
        feature = args.get("feature")
        top_n = args.get("top_n", default=20)

        # Ensure a gene name is provided:
        if not feature or feature == "":
            return Response(
                json.dumps({"error": "Please provide a feature(gene) name."}),
                status=400,
                mimetype="application/json",
            )

        h5_files_directory = (
            "compressed_data/h_sapiens/"  # Directory in the cloud storage bucket
        )

        # List all files in the bucket directory
        storage_client = storage.Client(project=os.getenv("GOOGLE_CLOUD_PROJECT"))
        files = storage_client.list_blobs(
            os.getenv("GOOGLE_CLOUD_BUCKET"), prefix=h5_files_directory
        )

        all_results = []

        # Loop through each .h5 file and get exp of the specified gene
        for file in files:
            if file.name.endswith(".h5"):
                try:
                    file_results = process_h5_file(
                        file.name, compute_gene_measurement, feature
                    )
                except ValueError as e:
                    # Return an error if the gene is not found in the gene mapping
                    return Response(
                        json.dumps({"error": str(e)}),
                        status=400,
                        mimetype="application/json",
                    )
                # Only extend all_results if file_results is not None
                if file_results is not None:
                    all_results.extend(file_results)
                else:
                    print(f"Skipping file {file.name} due to processing error.")

        # If there are no results, return an appropriate message
        if not all_results:
            return Response(
                json.dumps({"error": "No results found for the given gene."}),
                status=404,
                mimetype="application/json",
            )

        # Get the top N highest expressors from all the expression results
        result = get_highest_measurement(all_results, top_n)

        return Response(json.dumps(result), status=200, mimetype="application/json")
