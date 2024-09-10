import csv
import os
import json
import pandas as pd

from models.utils import convert_to_python_types


def load_gene_mapping():
    gene_csv_path = "data/ensembl_to_gene.csv"
    gene_to_ensembl = {}
    with open(gene_csv_path, newline="") as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header
        for row in reader:
            ensembl_id, gene_id = row
            gene_to_ensembl[gene_id] = ensembl_id
    return gene_to_ensembl


def compute_gene_measurement(adata, dataset_id, feature):
    """
    Compute all the gene expression measurement for a given feature (gene) from the AnnData object.

    Parameters:
        adata (AnnData): The AnnData object containing single-cell RNA-seq data.
        dataset_id (str): The identifier for the dataset (derived from file name).
        feature (str): The gene name (symbol) of interest.
        gene_mapping (dict): A dictionary mapping gene names to Ensembl IDs.

    Returns:
        list: A list of dictionaries containing information on cell type, tissue, disease, and expression levels of a gene.
    """

    ensembl_id = load_gene_mapping().get(feature.upper())
    if not ensembl_id:

        raise ValueError(f"Gene {feature} not found in our database.")
    
    if ensembl_id not in adata.var_names:
        return []  # Gene not found in this dataset

    gene_idx = adata.var_names.get_loc(ensembl_id)  # Get the index of the gene

    # Extract the expression data for the gene
    gene_expression = adata.layers["average"][:, gene_idx]
    fraction_detected = adata.layers["fraction"][:, gene_idx]
    
    try:
        with open(os.getenv("MANIFEST_FILE"), "r") as f:
            manifest = json.load(f)

        # Check if dataset_id exists in the manifest
        if dataset_id not in manifest:
            raise KeyError(f"Dataset ID {dataset_id} not found in the manifest.")
        
        metadata = manifest[dataset_id]

    except (KeyError, FileNotFoundError, json.JSONDecodeError) as e:
        # Log the error and skip this dataset
        print(f"Error processing metadata for {dataset_id}: {e}")
        return None  # Skip this dataset if there are any issues with metadata
    
    # Prepare results
    results = []
    for i in range(adata.n_obs):
        results.append(
            {
                "dataset_id": dataset_id,
                "dataset_title": metadata["dataset_title"],
                "collection_name": metadata["collection_name"],
                "cell_type": adata.obs["cell_type"][i],
                "tissue": adata.obs["tissue"][i],
                "disease": adata.obs["disease"][i],
                "sex": adata.obs["sex"][i],
                "development_stage": adata.obs["development_stage"][i],
                "average_expression": gene_expression[i],
                "fraction_detected": fraction_detected[i],
                "unit": metadata["unit"],
            }
        )
    return results


def get_highest_measurement(expressions, top_n):
    """
    Summarize and filter the top N highest expressors of a gene from all exp results.

    Parameters:
        expression (list): A list of expression measurements (dictionaries) across multiple datasets.
        top_n (int): Number of top expressors to return.

    Returns:
        list: A list of dictionaries containing the top N disease & cell types combination with the highest gene expression.
    """

    df = pd.DataFrame(expressions)
    avg_exp_df = (
        df[["disease", "cell_type", "unit", "collection_name", "average_expression"]]
        .groupby(by=["disease", "cell_type", "unit", "collection_name"])
        .mean()
        .reset_index()
    )

    avg_frac_df = (
        df[["disease", "cell_type", "unit", "collection_name", "fraction_detected"]]
        .groupby(by=["disease", "cell_type", "unit", "collection_name"])
        .mean()
        .reset_index()
    )

    combined = pd.merge(avg_exp_df, avg_frac_df, on=["collection_name", "disease", "cell_type", "unit"])

    sorted_results = sorted(
        combined.to_dict("records"), key=lambda x: x["average_expression"], reverse=True
    )
    final_result = sorted_results[:top_n]

    return convert_to_python_types(final_result)
