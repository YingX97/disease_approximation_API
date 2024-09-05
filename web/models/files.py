import os
from dotenv import load_dotenv
from google.cloud import storage
import scquill
from inspect import signature

load_dotenv()


def download_blob(bucket_name, source_blob_name, destination_file_name):
    #  Downloads a blob from the bucket if it does not already exist locally.
    if not os.path.exists(destination_file_name):
        storage_client = storage.Client(project=os.getenv("GOOGLE_CLOUD_PROJECT"))
        bucket = storage_client.bucket(bucket_name)
        blob = bucket.blob(source_blob_name)
        blob.download_to_filename(destination_file_name)


def process_h5_file(file_path, compute_func, *args):
    bucket_name = os.getenv("GOOGLE_CLOUD_BUCKET")
    local_path = f"/tmp/{os.path.basename(file_path)}"
    download_blob(bucket_name, file_path, local_path)

    app = scquill.Approximation()
    app = app.read_h5(local_path)
    
    # order of the group by columns needs to stay this way
    adata = app.to_anndata(
        groupby=(
            'cell_type', 'tissue', 'tissue_general', 
            'disease', 'sex', 'development_stage'
        )
    )

    dataset_id = os.path.basename(file_path).replace(".h5", "")

    # TO FIX: hard coded here, might not work in the future.
    # Check the number of arguments that the compute function expects
    sig = signature(compute_func)
    param_count = len(sig.parameters)

    return compute_func(adata, dataset_id, *args)
