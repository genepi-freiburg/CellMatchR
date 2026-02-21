import os
import numpy as np
from typing import Tuple
import pandas as pd
import logging
from huggingface_hub import hf_hub_download

from config import REQUIRED_FILES, DATAFOLDER, DATA_REPO_ID

logger = logging.getLogger(__name__)


def ensure_required_files():
    os.makedirs(DATAFOLDER, exist_ok=True)

    missing = [
        file for file in REQUIRED_FILES
        if not os.path.exists(os.path.join(DATAFOLDER, file))
    ]

    if missing:
        logger.warning(f"{len(missing)} Required file(s) are missing. Downloading them now...")
            
        download_huggingface_dataset(missing)


def load_test_data_settings(user_reference_datasets=None,
                            user_csv_path=None) -> dict:
    ensure_required_files()
    reference_data = []
    test_data = {}
    genelist = None

    #  Select reference datasets based on CLI argument
    required_files = REQUIRED_FILES

    if user_reference_datasets is not None:
        available_refs = {file.removesuffix(".parquet") for file, (datatype, _) in REQUIRED_FILES.items() if datatype == "reference"}
        unknown_refs = set(user_reference_datasets) - available_refs
        if unknown_refs:
            logger.error(f"Unknown reference dataset(s): {', '.join(sorted(unknown_refs))}")
            logger.info(f"Available: {', '.join(sorted(available_refs))}")
            exit(1)

        required_files = {
            file: (datatype, target_col) 
            for file, (datatype, target_col) in REQUIRED_FILES.items() 
            if not (datatype == "reference" and file.removesuffix(".parquet") not in user_reference_datasets)
        }

    logger.info(f"Using reference datasets: {', '.join(sorted(set(file.removesuffix('.parquet') for file, (datatype, _) in required_files.items() if datatype == 'reference')))}")

    for file, (datatype, target_col) in required_files.items():
        if datatype == "reference":
            parquet = pd.read_parquet(os.path.join(DATAFOLDER, file))
            reference_data.append(parquet)
        elif datatype == "test" and user_csv_path is None:
            df = pd.read_parquet(os.path.join(DATAFOLDER, file))
            test_data[file] = (df, target_col)
        elif datatype == "genelist":
            genelist = pd.read_csv(os.path.join(DATAFOLDER, file), header=None)[0]

    if user_csv_path is not None:
        csv_df = pd.read_csv(user_csv_path)
        csv_df.columns = csv_df.columns.str.upper()
        csv_name = os.path.splitext(os.path.basename(user_csv_path))[0]
        test_data[csv_name] = (csv_df, None)
        logger.info(f"Loaded user CSV '{csv_name}' with {csv_df.shape[0]} cells and {csv_df.shape[1]} genes.")

    reference_data = pd.concat(reference_data, axis=0, ignore_index=True)
    metadata_columns = [meta for meta in reference_data.columns if meta.startswith("meta_")]
    reference_genes = set(reference_data.columns) - set(metadata_columns)

    datasets = {}
    # For each test dataset, find the intersection of genes with the reference and genelist, and keep only those columns
    for name, (test_df, target_col) in test_data.items():
        test_genes = set(test_df.columns) - set(metadata_columns)
        gene_selection = sorted(reference_genes & test_genes & set(genelist))
        keep_cols = metadata_columns + gene_selection
        reference_data_subset = reference_data[keep_cols].copy()
        test_df_subset = test_df[[col for col in keep_cols if col in test_df.columns]].copy()
        if target_col is not None:
            reference_data_subset = reference_data_subset.rename(columns={target_col: "meta_target"})
            test_df_subset = test_df_subset.rename(columns={target_col: "meta_target"})
        else:
            reference_data_subset = reference_data_subset.rename(columns={"meta_celltype_tubular": "meta_target"})
            test_df_subset["meta_target"] = float("nan")
        datasets[name.removesuffix(".parquet")] = (reference_data_subset, test_df_subset)

    return datasets

def prepare_X_y(reference_data, test_data) -> Tuple[pd.DataFrame, pd.Series, pd.DataFrame, pd.Series | None]:
    meta_columns = [col for col in reference_data.columns if col.startswith("meta_")]
    X_train = reference_data.drop(columns=meta_columns)
    y_train = reference_data["meta_target"]
    X_test = test_data.drop(columns=meta_columns, errors="ignore")
    y_test = test_data["meta_target"]

    # Drop rows with missing target labels
    train_mask = ~y_train.isna()
    X_train = X_train[train_mask]
    y_train = y_train[train_mask]
    if y_test.isna().all():
        y_test = None
    else:
        test_mask = ~y_test.isna()
        X_test = X_test[test_mask]
        y_test = y_test[test_mask]

    # Log-transform the data
    X_train = np.log1p(X_train)
    X_test = np.log1p(X_test)
    return X_train, y_train, X_test, y_test

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