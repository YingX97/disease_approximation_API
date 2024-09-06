import json
import os
from dotenv import load_dotenv
from google.cloud import storage
from flask import Response, request
from flask_restful import Resource, abort

# helper functions that does the data processing
from models.highest_measurement import (
    compute_highest_measurement,
    zoom_in
)

from models.files import (
    process_h5_file
)

from api.v1.exceptions import (
    model_exceptions
)

load_dotenv()

class HighestMeasurement(Resource):
    
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
                mimetype="application/json"
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
                file_results = process_h5_file(
                    file.name,
                    compute_highest_measurement,
                    feature
                )
                all_results.extend(file_results)
        
        # Get the top N highest expressors from all the exp result
        result = zoom_in(all_results, top_n)
        
        return result