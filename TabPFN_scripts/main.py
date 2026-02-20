from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from tabpfn import TabPFNClassifier
import logging
from time import time
from torch.cuda import is_available
import numpy as np
from typing import Tuple

from sklearn.preprocessing import LabelEncoder


from utils import load_test_data_settings

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel("INFO")


def prepare_X_y(reference_data, test_data, target_col) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    meta_columns = [col for col in reference_data.columns if col.startswith("meta_")]
    X_train = reference_data.drop(columns=meta_columns)
    y_train = reference_data[target_col]
    X_test = test_data.drop(columns=meta_columns)
    y_test = test_data[target_col]

    # Drop rows with missing target labels
    train_mask = ~y_train.isna()
    X_train = X_train[train_mask]
    y_train = y_train[train_mask]
    test_mask = ~y_test.isna()
    X_test = X_test[test_mask]
    y_test = y_test[test_mask]

    # Log-transform the data
    X_train = X_train.apply(lambda x: np.log1p(x))
    X_test = X_test.apply(lambda x: np.log1p(x))
    return X_train, y_train, X_test, y_test

def fit_predict_evaluate(X_train, y_train, X_test, y_test):
    # Train TabPFN
    logger.info("Training TabPFN...")
    tabpfn_clf = TabPFNClassifier(device="cuda" if is_available() else "cpu")
    start_time = time()
    tabpfn_clf.fit(X_train.values, y_train.values)
    logger.info(f"TabPFN training time: {time() - start_time:.2f} seconds")

    # Train XGBoost
    logger.info("Training XGBoost...")
    # Encode labels for XGBoost
    le = LabelEncoder()
    y_train_encoded = le.fit_transform(y_train)
    y_test_encoded = le.transform(y_test)
    xgb_clf = XGBClassifier()
    start_time = time()
    xgb_clf.fit(X_train, y_train_encoded)
    logger.info(f"XGBoost training time: {time() - start_time:.2f} seconds")

    # Train Random Forest
    logger.info("Training Random Forest...")
    rf_clf = RandomForestClassifier()
    start_time = time()
    rf_clf.fit(X_train, y_train)
    logger.info(f"Random Forest training time: {time() - start_time:.2f} seconds")

    # Evaluate on test set
    tabpfn_preds = tabpfn_clf.predict(X_test.values)
    xgb_preds = xgb_clf.predict(X_test)
    rf_preds = rf_clf.predict(X_test)

    # Calculate and print accuracies
    tabpfn_acc = (tabpfn_preds == y_test.values).mean()
    xgb_acc = (xgb_preds == y_test_encoded).mean()
    rf_acc = (rf_preds == y_test.values).mean()

    logger.info(f"TabPFN Accuracy: {tabpfn_acc:.4f}")
    logger.info(f"XGBoost Accuracy: {xgb_acc:.4f}")
    logger.info(f"Random Forest Accuracy: {rf_acc:.4f}")

    return {"tabpfn": tabpfn_acc, "xgb": xgb_acc, "rf": rf_acc, "n_samples": len(y_test)}

def main():
    logger.info("Loading data...")
    test_settings = load_test_data_settings()

    target_col ="meta_celltype_tubular"

    for name, (reference_data, test_data) in test_settings.items():
        logger.info(f"Running on test dataset: {name}")

        X_train, y_train, X_test, y_test = prepare_X_y(reference_data, test_data, target_col)

        results = fit_predict_evaluate(X_train, y_train, X_test, y_test)

    for model in ["tabpfn", "xgb", "rf"]:
        total_correct = sum(r[model] * r["n_samples"] for r in results)
        total_samples = sum(r["n_samples"] for r in results)
        weighted_avg = total_correct / total_samples
        logger.info(f"{model}: {weighted_avg:.4f}")

if __name__ == "__main__":
    main()