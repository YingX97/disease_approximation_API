from collections import OrderedDict
import json
import os

from models.utils import (
    convert_to_python_types
)


def compute_diff_cell_abundance(adata, dataset_id, filters):
    """
    Computes the differential cell type abundance for a given dataset.

    Parameters:
        adata (AnnData): data obtained from scquill
        disease_keyword (str): Keyword to filter diseases.
        dataset_id (str): Identifier for the dataset.

    Returns:
        result (list): List of dictionaries containing the computed results.
    """
    filtered_obs = adata.obs
    
    # first check when disease keyword is given, make sure that the disease is present in adata
    # when disease is not present, we want to make sure that there is something not "normal" in adata
    if 'disease' not in filters:
        disease_filtered = filtered_obs[
            ~filtered_obs["disease"].str.contains('normal', case=False)
        ]["disease"].unique()
    else:
        disease_filtered = filtered_obs[
            filtered_obs["disease"].str.contains(str(filters['disease']), case=False)
        ]["disease"].unique()
    
    if len(disease_filtered) == 0:
        return []   
    
    # If a keyword is given, filter further
    for key in filters:
        if key not in ['unique_ids'] and filters[key] != '':
            if key == 'disease':
                filtered_obs = filtered_obs[
                    filtered_obs["disease"].str.contains(filters['disease'], case=False)
                    | (filtered_obs["disease"] == "normal")
                ]
            else:
                filtered_obs = filtered_obs[
                    filtered_obs[key].str.contains(str(filters[key]), case=False)
                ]
    
    if filtered_obs.empty:
        return []
    
    disease_status = filtered_obs["disease"].unique()
    cell_types = filtered_obs["cell_type"].unique()
    
    result = []
    with open(os.getenv("MANIFEST_FILE"), "r") as f:
        manifest = json.load(f)

    dataset_title = manifest[dataset_id]["dataset_title"]
    normal_den = filtered_obs[(filtered_obs["disease"] == "normal")]["cell_count"].sum()
    
    for cell_type in cell_types:
        for status in disease_status:
            if status == "normal":
                continue
            
            disease_den = filtered_obs[(filtered_obs["disease"] == status)]["cell_count"].sum()
            
            # disease_fraction = number of T cell cell / number of all cells (under the disease condition)
            normal_count = filtered_obs[
                (filtered_obs["cell_type"] == cell_type)
                & (filtered_obs["disease"] == "normal")
            ]["cell_count"].sum()
            disease_count = filtered_obs[
                (filtered_obs["cell_type"] == cell_type)
                & (filtered_obs["disease"] == status)
            ]["cell_count"].sum()
            
            if normal_den == 0 or disease_den == 0:
                continue
            
            normal_fraction = 1.0 * normal_count / normal_den
            disease_fraction = 1.0 * disease_count / disease_den
            delta_fraction = disease_fraction - normal_fraction

            # OrderedDict is used to ensure the dictionary is shown in the same order as defined
            result.append(
                OrderedDict(
                    [
                        ("disease", status),
                        ("dataset_title", dataset_title),
                        ("cell_type", cell_type),
                        ("comparison", "disease vs. normal"),
                        ("condition", "disease"),
                        ("condition_baseline", "normal"),
                        ("normal_count", normal_count),
                        ("disease_count", disease_count),
                        ("normal_total_count", normal_den),
                        ("disease_total_count", disease_den),
                        ("normal_fraction", normal_fraction),
                        ("disease_fraction", disease_fraction),
                        ("delta_fraction", delta_fraction),
                    ]
                )
            )

    return convert_to_python_types(result)