from collections import OrderedDict
import json
import os

import numpy as np

from models.utils import convert_to_python_types, load_ensembl_gene_pairs


def compute_diff_expression(
    adata, dataset_id, unit, log_scaled, N, filters
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
        adata.obs["disease"].str.contains(filters['disease'], case=False)
        | (adata.obs["disease"] == "normal")
    ]

    # If a keyword is given, filter further
    for key in filters:
        if key not in ['disease', 'unique_ids'] and filters[key] != '':
            filtered_obs = filtered_obs[
                filtered_obs[key].str.contains(str(filters[key]), case=False)
            ]

    if filtered_obs.empty:
        return []

    with open(os.getenv("MANIFEST_FILE"), "r") as f:
        manifest = json.load(f)

    dataset_title = manifest[dataset_id]["dataset_title"]
    disease_name = filtered_obs[
        filtered_obs["disease"].str.contains(filters['disease'], case=False)
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

            #cell type missing in either normal or disease condition will not be eligible for differential expression
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
