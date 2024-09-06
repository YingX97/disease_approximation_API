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

from models.differential_gene_expression import (
    compute_diff_expression
)

load_dotenv()

class DifferentialGeneExpression(Resource):
    """Get differential gene expression base on user filter"""
     
    @model_exceptions
    def post(self):
        disease_keyword = request.args.get("disease_keyword", default="", type=str)
        unique_ids = request.args.get("unique_ids", default="", type=str)  # comma-separated unique ids
        cell_type_keyword = request.args.get("cell_type", default="", type=str)
        sex_keyword = request.args.get("sex", default="", type=str)
        tissue_keyword = request.args.get("tissue", default="", type=str)
        stage_keyword = request.args.get("development_stage", default="", type=str)
        top_n = int(request.args.get("top_n", default=10, type=int))
        
        # Gather filters
        filters = {
            "disease": disease_keyword,
            "unique_ids": unique_ids.split(",") if unique_ids != '' else [],
            "cell_type": cell_type_keyword,
            "sex": sex_keyword,
            "tissue": tissue_keyword,
            "development_stage": stage_keyword
        }
        
        # Ensure only one of the keywords is provided, or neither
        if filters["disease"] != '' and len(filters['unique_ids']) > 0:
            return Response(
                json.dumps({"error": "Please provide either 'disease_keyword' or 'unique_ids', not both."}),
                status=400,
                mimetype="application/json"
            )
            
        # Ensure either disease_keyword or unique_ids is provided
        if filters["disease"] == '' and len(filters['unique_ids']) == 0:
            return Response(
                json.dumps({"error": "Either 'disease_keyword' or 'unique_ids' must be provided."}),
                status=400,
                mimetype="application/json"
            )

        matching_datasets = get_metadata(
            filters["disease"], 
            filters['cell_type'], 
            filters['unique_ids']
        )

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
                            compute_diff_expression,
                            d["unit"],
                            d["log_transformed"],
                            top_n,
                            filters
                        )
                        if result is not None:
                            all_results.extend(result)

        # return jsonify(all_results)
        response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
        return Response(response_json, mimetype="application/json")
