import os
import pandas as pd
import logging
from huggingface_hub import hf_hub_download

from config import REQUIRED_FILES, DATAFOLDER, DATA_REPO_ID

logger = logging.getLogger(__name__)


def ensure_required_files():
    if not os.path.exists(DATAFOLDER):
        os.makedirs(DATAFOLDER)

    missing = []
    for file in REQUIRED_FILES:
        if not os.path.exists(os.path.join(DATAFOLDER, file)):
            missing.append(file)
    
    if len(missing) > 0:
        logger.warning(f"{len(missing)} Required file(s) are missing. Downloading them now...")
            
        download_huggingface_dataset(missing)


def load_test_data_settings() -> dict:
    ensure_required_files()
    reference_data = []
    test_data = {}
    genelist = None

    for file, datatype in REQUIRED_FILES.items():
        if datatype == "reference":
            parquet = pd.read_parquet(os.path.join(DATAFOLDER, file))
            reference_data.append(parquet)
        elif datatype == "test":
            test_data[file] = pd.read_parquet(os.path.join(DATAFOLDER, file))
        elif datatype == "genelist":
            genelist = pd.read_csv(os.path.join(DATAFOLDER, file), header=None)[0]

    reference_data = pd.concat(reference_data, axis=0, ignore_index=True)
    metadata_columns = [meta for meta in reference_data.columns if meta.startswith("meta_")]
    reference_genes = set(reference_data.columns) - set(metadata_columns)

    datasets = {}
    # For each test dataset, find the intersection of genes with the reference and genelist, and keep only those columns
    for name, test_df in test_data.items():
        test_genes = set(test_df.columns) - set(metadata_columns)
        gene_selection = sorted(reference_genes & test_genes & set(genelist))
        keep_cols = metadata_columns + gene_selection
        datasets[name] = (reference_data[keep_cols], test_df[keep_cols])

    return datasets


def download_huggingface_dataset(missing):
    
    logger.info(f"Downloading {len(missing)} missing file(s)...")
    for filename in missing:
        hf_hub_download(
            repo_id=DATA_REPO_ID,
            filename=filename,
            local_dir=DATAFOLDER,
            repo_type="dataset",
        )
    logger.info("Download complete.")