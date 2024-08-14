import os
import re
import json
import numpy as np
from dotenv import load_dotenv
from collections import OrderedDict
from utils.helpers import convert_to_python_types, load_ensembl_gene_pairs

load_dotenv()


def compute_diff_cell_abundance(adata, disease_keyword, dataset_id):
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
        df_obs["disease"].str.contains(disease_keyword, case=False)
        | (df_obs["disease"] == "normal")
    ]
    disease_name = filtered_obs[
        filtered_obs["disease"].str.contains(disease_keyword, case=False)
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


def compute_diff_expression(
    adata, disease_keyword, dataset_id, unit, log_scaled, N, cell_type_keyword
):
    """
    Computes the differential cell type abundance for a given dataset.

    Parameters:
        adata (AnnData): data obtained from scquill
        keyword (str): Keyword to filter diseases.
        dataset_id (str): Identifier for the dataset.

    Returns:
        result (list): List of dictionaries containing the computed results.
    """
    # Ensure ENSEMBL_TO_GENE_FILE is loaded correctly
    ensembl_gene_pair_file = os.getenv("ENSEMBL_TO_GENE_FILE")

    # Load Ensembl to gene name mapping
    ensembl_to_gene = load_ensembl_gene_pairs(ensembl_gene_pair_file)

    expression_data = adata.layers["average"]
    fraction_data = adata.layers["fraction"]

    # Set numerical indices on the original adata.obs
    adata.obs["numerical_index"] = np.arange(adata.obs.shape[0])

    # Filter the obs to only include rows with the disease_keyword or 'normal'
    filtered_obs = adata.obs[
        adata.obs["disease"].str.contains(disease_keyword, case=False)
        | (adata.obs["disease"] == "normal")
    ]

    # If a cell type keyword is given, filter further

    if cell_type_keyword:
        filtered_obs = filtered_obs[
            filtered_obs["cell_type"].str.contains(cell_type_keyword, case=False)
        ]

    if filtered_obs.empty:
        return []

    with open(os.getenv("MANIFEST_FILE"), "r") as f:
        manifest = json.load(f)

    dataset_title = manifest[dataset_id]["dataset_title"]
    disease_name = filtered_obs[
        filtered_obs["disease"].str.contains(disease_keyword, case=False)
    ]["disease"].unique()[0]

    # Convert filtered_obs['numerical_index'] to a numpy array of integer indices
    integer_indices = filtered_obs["numerical_index"].to_numpy()

    # Filter the expression and fraction data using the filtered observations' integer indices
    filtered_expression_data = expression_data[integer_indices, :]
    filtered_fraction_data = fraction_data[integer_indices, :]

    cell_types = filtered_obs["cell_type"].unique()
    disease_status = filtered_obs["disease"].unique()
    result = []

    for cell_type in cell_types:

        for status in disease_status:
            if status == "normal":
                continue

            # Get the numerical indices for normal and disease
            cell_type_in_normal = filtered_obs[
                (filtered_obs["cell_type"] == cell_type)
                & (filtered_obs["disease"] == "normal")
            ]["numerical_index"]
            cell_type_in_disease = filtered_obs[
                (filtered_obs["cell_type"] == cell_type)
                & (filtered_obs["disease"] == status)
            ]["numerical_index"]

            # cell type missing in either normal or disease condition will not be eligible for differential expression
            if len(cell_type_in_normal) == 0 or len(cell_type_in_disease) == 0:
                continue

            normal_idx = cell_type_in_normal.values[0]
            disease_idx = cell_type_in_disease.values[0]

            # Map the original indices to the filtered data indices
            normal_idx = np.where(integer_indices == normal_idx)[0][0]
            disease_idx = np.where(integer_indices == disease_idx)[0][0]

            # Ensure normal_idx and disease_idx are within bounds
            if normal_idx >= len(filtered_expression_data) or disease_idx >= len(
                filtered_expression_data
            ):
                continue

            normal_expr = filtered_expression_data[normal_idx, :]
            disease_expr = filtered_expression_data[disease_idx, :]
            normal_fraction = filtered_fraction_data[normal_idx, :]
            disease_fraction = filtered_fraction_data[disease_idx, :]
            delta_fraction = disease_fraction - normal_fraction

            log2_fc = np.log2((disease_expr + 1) / (normal_expr + 1))

            # ensure top_n does not exceed the number of features
            if N >= log2_fc.shape[0]:
                N = log2_fc.shape[0] - 1

            # Sort by delta fraction
            top_up_indices = np.argsort(delta_fraction)[-N:][::-1]
            top_down_indices = np.argsort(delta_fraction)[:N]

            # Combine up and down indices with regulation labels
            combined_indices = [(idx, "up") for idx in top_up_indices] + [
                (idx, "down") for idx in top_down_indices
            ]

            for idx, regulation in combined_indices:
                gene_name = ensembl_to_gene.get(
                    adata.var_names[idx], adata.var_names[idx]
                )
                result.append(
                    OrderedDict(
                        [
                            ("disease", disease_name),
                            ("dataset_title", dataset_title),
                            ("cell_type", cell_type),
                            ("comparison", "disease vs. normal"),
                            ("condition", "disease"),
                            ("condition_baseline", "normal"),
                            ("regulation", regulation),
                            ("gene", gene_name),
                            ("unit", unit),
                            ("normal_expr", normal_expr[idx]),
                            ("disease_expr", disease_expr[idx]),
                            ("log_scaled", log_scaled),
                            ("log2_fc", log2_fc[idx]),
                            ("normal_fraction", normal_fraction[idx]),
                            ("disease_fraction", disease_fraction[idx]),
                            ("delta_fraction", delta_fraction[idx]),
                        ]
                    )
                )

    return convert_to_python_types(result)


def get_metadata(disease_keyword="", cell_type_keyword="", unique_ids=[]):
    result = []
    with open(os.getenv("MANIFEST_FILE"), "r") as f:
        manifest = json.load(f)

    # Filter by cell type keyword if provided
    if cell_type_keyword != "":
        pattern = re.compile(
            rf"\b{re.escape(cell_type_keyword.lower())}\b", re.IGNORECASE
        )
        for dataset_id in manifest:
            metadata = manifest[dataset_id]
            if "cell_types" in metadata:
                # Exact match, considering word boundaries
                matching_cell_types = [
                    cell_type
                    for cell_type in metadata["cell_types"]
                    if pattern.search(cell_type.lower())
                ]
                # Only append the item if there is a valid match with the cell_type_keyword
                if matching_cell_types:
                    item = {
                        "uid": metadata["ids"][0],
                        "disease": metadata["diseases"],
                        "dataset_id": dataset_id,
                        "dataset_title": metadata["dataset_title"],
                        "collection_name": metadata["collection_name"],
                        # "matched_cell_types": matching_cell_types,
                        "cell_types": metadata["cell_types"],
                        "unit": metadata["unit"],
                        "log_scaled": metadata["log_scaled"],
                        "has_normal_baseline": metadata["has_normal_baseline"],
                    }
                    result.append(item)
    # Filter by disease keyword if provided
    else:
        for dataset_id in manifest:
            metadata = manifest[dataset_id]
            if "diseases" in metadata:
                for unique_id, disease in zip(metadata["ids"], metadata["diseases"]):
                    item = {
                        "uid": unique_id,
                        "disease": disease,
                        "dataset_id": dataset_id,
                        "dataset_title": metadata["dataset_title"],
                        "collection_name": metadata["collection_name"],
                        "cell_types": metadata["cell_types"],
                        "unit": metadata["unit"],
                        "log_scaled": metadata["log_scaled"],
                        "has_normal_baseline": metadata["has_normal_baseline"],
                    }
                    if len(unique_ids) > 0 and unique_ids[0] != "":
                        if unique_id in unique_ids:
                            result.append(item)
                    else:
                        if disease_keyword.lower() in disease.lower():
                            result.append(item)

    return convert_to_python_types(result)
