from flask import Blueprint, request, Response
from dotenv import load_dotenv
from google.cloud import storage
from utils.file_handling import process_h5_file
from utils.data_preprocessing import get_metadata
import json
import os

load_dotenv()

get_metadata_bp = Blueprint("metadata", __name__)


@get_metadata_bp.route("/metadata", methods=["GET"])
def get_metadata_route():
    disease_keyword = request.args.get("disease_keyword", default="", type=str)
    cell_type_keyword = request.args.get("cell_type_keyword", default="", type=str)

    # Ensure only one of the keywords is provided, or neither
    if disease_keyword and cell_type_keyword:
        return Response(
            json.dumps({"error": "Please provide either 'disease_keyword' or 'cell_type_keyword', not both."}),
            status=400,
            mimetype="application/json"
        )

    # Ensure at least one keyword is provided
    if not disease_keyword and not cell_type_keyword:
        return Response(
            json.dumps({"error": "Please provide either 'disease_keyword' or 'cell_type_keyword'."}),
            status=400,
            mimetype="application/json"
        )
    all_results = get_metadata(disease_keyword, cell_type_keyword)

    # Convert the list of OrderedDict to JSON string to preserve the order
    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype="application/json")
