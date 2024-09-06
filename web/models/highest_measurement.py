import csv
import os
import pandas as pd

from models.utils import (
    convert_to_python_types
)

def load_gene_mapping():
    gene_csv_path = "data/ensembl_to_gene.csv"
    gene_to_ensembl = {}
    with open(gene_csv_path, newline='') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)  # Skip the header
        for row in reader:
            ensembl_id, gene_id = row
            gene_to_ensembl[gene_id] = ensembl_id
    return gene_to_ensembl

def compute_highest_measurement(adata, dataset_id, feature):
    
    ensembl_id = load_gene_mapping().get(feature.upper())
    if not ensembl_id:
    
        raise ValueError(f"Gene {feature} not found in the gene mapping.")
    if ensembl_id not in adata.var_names:
        return []  # Gene not found in this dataset
    
    gene_idx = adata.var_names.get_loc(ensembl_id)  # Get the index of the gene

    # Extract the expression data for the gene
    gene_expression = adata.layers['average'][:, gene_idx]
    fraction_detected = adata.layers['fraction'][:, gene_idx]

    # Prepare results
    results = []
    for i in range(adata.n_obs):
        results.append({
            'dataset_id': dataset_id,
            'cell_type': adata.obs['cell_type'][i],
            'tissue': adata.obs['tissue'][i],
            'disease': adata.obs['disease'][i],
            'sex': adata.obs['sex'][i],
            'development_stage': adata.obs['development_stage'][i],
            'average_expression': gene_expression[i],
            'fraction_detected': fraction_detected[i]
        })
    
    return results

def zoom_in(expressions, top_n):
    df = pd.DataFrame(expressions)
    sum_avg_df = df[['disease', 'cell_type', 'average_expression']]\
        .groupby(by=['disease', 'cell_type'])\
        .sum().reset_index()  
    
    avg_frac_df = df[['disease', 'cell_type', 'fraction_detected']]\
        .groupby(by=['disease', 'cell_type'])\
        .mean().reset_index()  
    
    combined = pd.merge(sum_avg_df, avg_frac_df, on=['disease', 'cell_type'])
    
    sorted_results = sorted(combined.to_dict('records'), key=lambda x: x['average_expression'], reverse=True)
    final_result = sorted_results[:top_n]
    
    return convert_to_python_types(final_result)