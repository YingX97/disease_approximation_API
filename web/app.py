"""
Web application supporting the cell atlas disease approximation API
"""
from flask import (
    Flask,
)
from flask_restful import Api
from flask_cors import CORS
from config import configuration as config
from api import api_dict


##############################
app = Flask(__name__, static_url_path="/static", template_folder="templates")
app_api = Api(app)
with open('secret_key.txt') as f:
    app.config['SECRET_KEY'] = f.read()
##############################

api_authorized_versions = [config['api_version']]
authorized_resources = {}
for api_version in api_authorized_versions:
    # Connect to endpoints
    if config["api_version"] == "v1":
        get_api_endpoint = api_dict[api_version]["endpoint_handler"]
        for api_name, api_object in api_dict[api_version]["objects"].items():
            app_api.add_resource(api_object, get_api_endpoint(api_name))

    # Open to cross-origin requests
    authorized_resources[f"/{api_version}/*"] = {"origins": "*"}

# Cross-origin request handler
CORS(app, resources=authorized_resources)


# Main loop
if __name__ == "__main__":
    # app.run(debug=True)
    app.run(host="0.0.0.0", port=5000, debug=True)