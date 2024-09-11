import os
import re
import json
from models.utils import *


def get_metadata(disease_keyword="", cell_type_keyword="", unique_ids=[]):
    result = []

    with open(os.getenv("MANIFEST_FILE"), "r") as f:
        manifest = json.load(f)

    # if unique id is provided, ignore disease and cell type keyword
    if len(unique_ids) > 0 and unique_ids[0] != "":
        for dataset_id in manifest:
            metadata = manifest[dataset_id]
            for unique_id in metadata["ids"]:
                if unique_id in unique_ids:
                    item = {
                        "uid": metadata["ids"],
                        "disease": metadata["disease"],
                        "dataset_id": dataset_id,
                        "dataset_title": metadata["dataset_title"],
                        "collection_name": metadata["collection_name"],
                        "cell_types": metadata["cell_type"],
                        "unit": metadata["unit"],
                        "log_transformed": metadata["log_transformed"],
                        "has_normal_baseline": metadata["has_normal_baseline"],
                    }
                    result.append(item)

        return convert_to_python_types(result)

    # if user specified a cell type keyword, e.g "NK"
    if cell_type_keyword != "" and disease_keyword == "":
        pattern = re.compile(
            rf"\b{re.escape(cell_type_keyword.lower())}\b", re.IGNORECASE
        )
        for dataset_id in manifest:
            metadata = manifest[dataset_id]
            if "cell_type" in metadata:
                matching_cell_types = [
                    cell_type
                    for cell_type in metadata["cell_type"]
                    if pattern.search(cell_type.lower())
                ]
                # Only append the item if there is a valid match with the cell_type_keyword
                # associate each unique_id with it's disease
                if len(matching_cell_types) > 0:
                    for i, uid in enumerate(metadata["ids"]):
                        disease = metadata["disease"][i]
                        item = {
                            "uid": uid,
                            "disease": disease,
                            "dataset_id": dataset_id,
                            "dataset_title": metadata["dataset_title"],
                            "collection_name": metadata["collection_name"],
                            "cell_types": metadata["cell_type"],
                            "unit": metadata["unit"],
                            "log_transformed": metadata["log_transformed"],
                            "has_normal_baseline": metadata["has_normal_baseline"],
                        }
                        result.append(item)

    # if user specified a disease keyword, e.g "covid"
    elif disease_keyword != "" and cell_type_keyword == "":
        for dataset_id in manifest:
            metadata = manifest[dataset_id]
            if "disease" in metadata:
                disease_match = [
                    d
                    for d in metadata["disease"]
                    if disease_keyword.lower() in d.lower()
                ]
                if len(disease_match) > 0:
                    for i, disease in enumerate(disease_match):
                        uid = metadata["ids"][i],
                        item = {
                            "uid": uid,
                            "disease": disease,
                            "dataset_id": dataset_id,
                            "dataset_title": metadata["dataset_title"],
                            "collection_name": metadata["collection_name"],
                            "cell_types": metadata["cell_type"],
                            "unit": metadata["unit"],
                            "log_transformed": metadata["log_transformed"],
                            "has_normal_baseline": metadata["has_normal_baseline"],
                        }
                        result.append(item)

    # if user provide both disease and cell type, e.g "covid that contains CD8-positive T cell"
    elif disease_keyword != "" and cell_type_keyword != "":
        pattern = re.compile(
            rf"\b{re.escape(cell_type_keyword.lower())}\b", re.IGNORECASE
        )
        for dataset_id in manifest:
            metadata = manifest[dataset_id]
            if "cell_type" in metadata and "disease" in metadata:
                # Exact match, considering word boundaries
                matching_cell_types = [
                    cell_type
                    for cell_type in metadata["cell_type"]
                    if pattern.search(cell_type.lower())
                ]

                disease_match = [
                    d
                    for d in metadata["disease"]
                    if disease_keyword.lower() in d.lower()
                ]

                # Only append the item if there is a valid match with the cell_type_keyword
                if len(matching_cell_types) > 0 and len(disease_match) > 0:
                    for i, disease in enumerate(disease_match):
                        uid = metadata["ids"][i],
                        item = {
                            "uid": uid,
                            "disease": disease,
                            "dataset_id": dataset_id,
                            "dataset_title": metadata["dataset_title"],
                            "collection_name": metadata["collection_name"],
                            "cell_types": metadata["cell_type"],
                            "unit": metadata["unit"],
                            "log_transformed": metadata["log_transformed"],
                            "has_normal_baseline": metadata["has_normal_baseline"],
                        }
                        result.append(item)
    # if no keyword provided, show metadata for all available disease
    else:
        for dataset_id in manifest:
            metadata = manifest[dataset_id]
            for i, uid in enumerate(metadata["ids"]):
                disease = metadata["disease"][i] 
                item = {
                    "uid": uid,
                    "disease": disease,
                    "dataset_id": dataset_id,
                    "dataset_title": metadata["dataset_title"],
                    "collection_name": metadata["collection_name"],
                    "cell_types": metadata["cell_type"],
                    "tissues": metadata["tissue_general"],
                    "unit": metadata["unit"],
                    "log_transformed": metadata["log_transformed"],
                    "has_normal_baseline": metadata["has_normal_baseline"],
                }
                result.append(item)

    return convert_to_python_types(result)
