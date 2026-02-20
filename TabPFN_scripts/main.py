from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier
from tabpfn import TabPFNClassifier
import logging
from time import time
from torch.cuda import is_available

from sklearn.preprocessing import LabelEncoder


from utils import load_data


logger = logging.getLogger(__name__)
logger.setLevel("INFO")

def main():
    logger.info("Loading data...")
    reference_data, test_data = load_data()

    meta_columns = [col for col in reference_data.columns if col.startswith("meta_")]
    target_col ="meta_celltype_tubular"

    X_train = reference_data.drop(columns=meta_columns)
    y_train = reference_data[target_col]

    # Drop rows with missing target labels in training data
    train_mask = ~y_train.isna()
    X_train = X_train[train_mask]
    y_train = y_train[train_mask]

    logger.info(f"Training models on {X_train.shape[0]} samples and {X_train.shape[1]} features...")

    X_test = test_data.drop(columns=meta_columns)
    y_test = test_data[target_col]

    logger.info(f"Evaluating models on {X_test.shape[0]} samples and {X_test.shape[1]} features...")

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


if __name__ == "__main__":
    main()