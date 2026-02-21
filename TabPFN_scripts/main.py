
from tabpfn import TabPFNClassifier
import logging
from time import time
from torch.cuda import is_available
import numpy as np
from typing import Tuple
import argparse
import os 
import pandas as pd
import matplotlib.pyplot as plt


from utils import load_test_data_settings

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)
logger.setLevel("INFO")


def plot_probabilities(probs: pd.DataFrame, dataset_name: str):
    n_total = len(probs)
    n_plots = (n_total + 8) // 9

    for plot in range(n_plots):
        batch = probs.iloc[plot * 9 : (plot + 1) * 9]
        fig, axes = plt.subplots(3, 3, figsize=(12, 8))

        for ax, (idx, row) in zip(axes.flat, batch.iterrows()):
            row.plot.bar(ax=ax)
            ax.set_ylabel("Probability")
            ax.set_title(f"Sample {idx}")

        # Hide unused subplots
        for ax in axes.flat[len(batch):]:
            ax.set_visible(False)

        fig.tight_layout()
        fig.savefig(f"results/probabilities_{dataset_name}_part{plot + 1}.png")
        plt.close(fig)

        

def prepare_X_y(reference_data, test_data) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    meta_columns = [col for col in reference_data.columns if col.startswith("meta_")]
    X_train = reference_data.drop(columns=meta_columns)
    y_train = reference_data["meta_target"]
    X_test = test_data.drop(columns=meta_columns)
    y_test = test_data["meta_target"]

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
    tabpfn_clf = TabPFNClassifier(device="cuda" if is_available() else "cpu",
                                  ignore_pretraining_limits=True) # Ignoring limits for local CPU training
    start_time = time()
    tabpfn_clf.fit(X_train.values, y_train.values)
    logger.info(f"TabPFN training time: {time() - start_time:.2f} seconds")

    # Evaluate on test set
    logger.info("Predicting with TabPFN...")
    start_time = time()
    tabpfn_probs = tabpfn_clf.predict_proba(X_test.values)
    logger.info(f"TabPFN prediction time: {time() - start_time:.2f} seconds")

    tabpfn_probs = pd.DataFrame(tabpfn_probs, columns=tabpfn_clf.classes_, index=X_test.index)

    # Calculate and print accuracies
    tabpfn_preds = tabpfn_probs.idxmax(axis=1)

    if y_test is not None:
        tabpfn_acc = (tabpfn_preds == y_test).mean()

    logger.info(f"TabPFN Accuracy: {tabpfn_acc:.4f}")

    return {
        "acc": tabpfn_acc if y_test is not None else None,
        "n_samples": len(y_test) if y_test is not None else 0,
        "preds": tabpfn_preds,
        "probs": tabpfn_probs
        }

def main():
    # Setup CLI
    parser = argparse.ArgumentParser(description="Evaluate TabPFN test datasets.")

    parser.add_argument("--reference_datasets", nargs="+", default=None,
                        help="List of reference datasets to use (default: all)")
    
    args = parser.parse_args()


    logger.info("Loading data...")
    test_settings = load_test_data_settings(args.reference_datasets)

    os.makedirs("results", exist_ok=True)

    all_results = []
    for name, (reference_data, test_data) in test_settings.items():
        logger.info("=" * 60)
        logger.info(f"Test dataset: {name}")
        logger.info("-" * 60)

        X_train, y_train, X_test, y_test = prepare_X_y(reference_data, test_data)

        results = fit_predict_evaluate(X_train, y_train, X_test, y_test)

        plot_probabilities(results["probs"], name)
        logger.info(f"Saved probability plot for {name}.")

        pd.DataFrame(results["probs"]).to_csv(f"results/results_{name}.csv", index=False)
        logger.info(f"Saved probabilities for {name} to CSV.")

        all_results.append(results)

    # This can only be calcualted it we have true labels for the test set
    if all_results[0]["acc"] is not None:
        logger.info("=" * 60)
        logger.info("WEIGHTED AVERAGE ACCURACY")
        logger.info("-" * 60)
        total_correct = sum(r["acc"] * r["n_samples"] for r in all_results)
        total_samples = sum(r["n_samples"] for r in all_results)
        weighted_avg = total_correct / total_samples if total_samples > 0 else 0
        logger.info(f"TabPFN: {weighted_avg:.4f}")

   

if __name__ == "__main__":
    main()