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
    df_obs = adata.obs
    filtered_obs = df_obs[
        df_obs["disease"].str.contains(filters['disease'], case=False)
        | (df_obs["disease"] == "normal")
    ]
    
    # If a keyword is given, filter further
    for key in filters:
        if key not in ['disease', 'unique_ids'] and filters[key] != '':
            filtered_obs = filtered_obs[
                filtered_obs[key].str.contains(str(filters[key]), case=False)
            ]

    if filtered_obs.empty:
        return []
    
    disease_name = filtered_obs[
        filtered_obs["disease"].str.contains(filters['disease'], case=False)
    ]["disease"].unique()[0]
    result = []
    with open(os.getenv("MANIFEST_FILE"), "r") as f:
        manifest = json.load(f)

    dataset_title = manifest[dataset_id]["dataset_title"]
    normal_den = filtered_obs[(filtered_obs["disease"] == "normal")]["cell_count"].sum()
    disease_den = filtered_obs[(filtered_obs["disease"] != "normal")][
        "cell_count"
    ].sum()
    for cell_type in filtered_obs.cell_type.unique():
        # disease_fraction = number of T cell cell / number of all cells (under the disease condition)
        normal_count = filtered_obs[
            (filtered_obs["cell_type"] == cell_type)
            & (filtered_obs["disease"] == "normal")
        ]["cell_count"].sum()
        disease_count = filtered_obs[
            (filtered_obs["cell_type"] == cell_type)
            & (filtered_obs["disease"] != "normal")
        ]["cell_count"].sum()
        normal_fraction = 1.0 * normal_count / normal_den
        disease_fraction = 1.0 * disease_count / disease_den
        delta_fraction = disease_fraction - normal_fraction

        # OrderedDict is used to ensure the dictionary is shown in the same order as defined
        result.append(
            OrderedDict(
                [
                    ("disease", disease_name),
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