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
    """
        Get a list of metadata based on disease, cell type, or tissue keywords
        This quickly help answer question such as:
        - is  there any record of the disease I am interested in?
        - What disease and datasets will contains the my cell type of interes?
        - What disease and datasets will contains the tissue I want?
    """
    
    @model_exceptions
    def get(self):
        disease_keyword = request.args.get("disease", default="", type=str)
        cell_type_keyword = request.args.get("cell_type", default="", type=str)
        tissue_keyword = request.args.get("tissue", default="", type=str)
        
        all_results = get_metadata(disease_keyword, cell_type_keyword, tissue_keyword)
        
        # Convert the list of OrderedDict to JSON string to preserve the order
        response_json = json.dumps(all_results, ensure_ascii=False, indent=4)
        return Response(response_json, mimetype="application/json")
            