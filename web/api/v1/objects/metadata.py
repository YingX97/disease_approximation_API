import json
from flask import Response, request
from flask_restful import Resource

from models.metadata import (
    get_metadata
)

from api.v1.exceptions import (
    model_exceptions
)

class Metadata(Resource):
    """Get a list of metadata based on user's../."""

    
    @model_exceptions
    def get(self):
        disease_keyword = request.args.get("disease_keyword", default="", type=str)
        cell_type_keyword = request.args.get("cell_type_keyword", default="", type=str)
        
        all_results = get_metadata(disease_keyword, cell_type_keyword)
        
        # Convert the list of OrderedDict to JSON string to preserve the order
        response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
        return Response(response_json, mimetype="application/json")
            