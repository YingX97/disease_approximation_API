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

load_dotenv()

class HighestExpressor(Resource):
    
    @model_exceptions
    def post(self):
        return {'message': 'Endpoint for highest expression'}