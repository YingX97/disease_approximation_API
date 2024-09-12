import json
import os
from dotenv import load_dotenv
from google.cloud import storage
from flask import Response, request
from flask_restful import Resource, abort

from models.metadata import (
    get_metadata
)

from models.files import (
    process_h5_file
)

from api.v1.exceptions import (
    model_exceptions
)

from models.differential_cell_type_abundance import (
    compute_diff_cell_abundance
)

class DifferentialCellTypeAbundance(Resource):
    """
    Get differential cell type abundance for a disease or tissue of interest across multiple datasets.
    
    This API helps answer questions such as:
    1. "In COVID-19, what is the differential cell abundance of each cell type between normal and disease states?"
    2. "Across all diseases affecting the lung, what is the differential cell abundance of each cell type in that tissue between normal and disease states?"
    """
    
    @model_exceptions
    def post(self):
        disease_keyword = request.args.get("disease", default="", type=str)
        cell_type_keyword = request.args.get("cell_type", default="", type=str)
        # Get unique_ids as a comma-separated string and split it into a list
        unique_ids_str = request.args.get("unique_ids", default="", type=str)
        unique_ids = unique_ids_str.split(",") if unique_ids_str else []
            
        # Gather filters
        filters = {
            "disease": disease_keyword,
            "cell_type": cell_type_keyword,
            "unique_ids": unique_ids,
        }

        matching_datasets = get_metadata(
            filters["disease"], 
            filters["cell_type"], 
            "",
            filters["unique_ids"]
        )

        (matching_datasets)
        if len(matching_datasets) == 0:
            return {"message": "No datasets found satisfying disease and cell type filter"}

        h5_files_directory = (
            "compressed_data/h_sapiens/"  # Directory in the cloud storage bucket
        )

        all_results = []

        # List all files in the bucket directory
        storage_client = storage.Client(project=os.getenv("GOOGLE_CLOUD_PROJECT"))
        files = storage_client.list_blobs(
            os.getenv("GOOGLE_CLOUD_BUCKET"), prefix=h5_files_directory
        )
        
        for file in files:
            if file.name.endswith(".h5"):
                dataset_id = str(file.name).split("/")[-1].replace(".h5", "")
                for index, d in enumerate(matching_datasets):
                    if dataset_id == d["dataset_id"]:
                        result = process_h5_file(
                            file.name,
                            compute_diff_cell_abundance,
                            filters
                        )
                        if result is not None:
                            all_results.extend(result)

        # Convert the list of OrderedDict to JSON string to preserve the order
        response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
        return Response(response_json, mimetype="application/json")
