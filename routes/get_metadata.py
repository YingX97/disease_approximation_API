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
    all_results = get_metadata(disease_keyword)

    # Convert the list of OrderedDict to JSON string to preserve the order
    response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
    return Response(response_json, mimetype="application/json")
