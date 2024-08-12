import numpy as np
import pandas as pd


def convert_to_python_types(result):
    """
    Converts NumPy types in a list of dictionaries to native Python types.

    Parameters:
        result (list): List of dictionaries containing results.

    Returns:
        list: List of dictionaries with values converted to Python types.

    This function fixes the error: TypeError: Object of type int64 is not JSON serializable.
    """
    for entry in result:
        for key, value in entry.items():
            if isinstance(value, (np.integer, np.int64)):
                entry[key] = int(value)
            elif isinstance(value, (np.floating, np.float64)):
                entry[key] = float(value)
            elif isinstance(value, (np.ndarray, np.generic)):
                entry[key] = value.tolist()
    return result


def load_ensembl_gene_pairs(filename):
    """
    Loads Ensembl to gene name mapping from a CSV file.

    Parameters:
        filename (str): Path to the CSV file.

    Returns:
        dict: A dictionary mapping Ensembl IDs to gene names.
    """
    df = pd.read_csv(filename)
    return pd.Series(df.gene_id.values, index=df.ensembl_id).to_dict()
