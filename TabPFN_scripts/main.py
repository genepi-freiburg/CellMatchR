import argparse
import logging
import os
from time import time
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tabpfn import TabPFNClassifier
from torch.cuda import is_available

from utils import load_test_data_settings, prepare_X_y

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s')
logger = logging.getLogger(__name__)


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

    tabpfn_acc = None
    if y_test is not None:
        tabpfn_acc = (tabpfn_preds == y_test).mean()
        logger.info(f"TabPFN Accuracy: {tabpfn_acc:.4f}")

    return {
        "acc": tabpfn_acc,
        "n_samples": len(y_test) if y_test is not None else 0,
        "preds": tabpfn_preds,
        "probs": tabpfn_probs
        }

def main():
    # Setup CLI
    parser = argparse.ArgumentParser(description="Evaluate TabPFN test datasets.")

    parser.add_argument("--reference_datasets", nargs="+", default=None,
                        help="List of reference datasets to use (default: all)")
    parser.add_argument("--csv", type=str, default=None,
                        help="Path to a CSV file with cells as rows and genes as columns (gene names as headers)")

    args = parser.parse_args()


    logger.info("Loading data...")
    test_settings = load_test_data_settings(args.reference_datasets, user_csv_path=args.csv)

    os.makedirs("results", exist_ok=True)

    all_results = []
    for name, (reference_data, test_data) in test_settings.items():
        logger.info("=" * 60)
        logger.info(f"Test dataset: {name}")
        logger.info("-" * 60)
        
        X_train, y_train, X_test, y_test = prepare_X_y(reference_data, test_data)
        logger.info(f"After gene selection {X_test.shape[1]} genes remain for matching.")

        results = fit_predict_evaluate(X_train, y_train, X_test, y_test)

        plot_probabilities(results["probs"], name)
        logger.info(f"Saved probability plot for {name}.")

        results["probs"].to_csv(f"results/results_{name}.csv", index=False)
        logger.info(f"Saved probabilities for {name} to CSV.")

        all_results.append(results)

    # This can only be calculated if we have true labels for the test set
    if any(r["acc"] is not None for r in all_results):
        logger.info("=" * 60)
        logger.info("WEIGHTED AVERAGE ACCURACY")
        logger.info("-" * 60)
        total_correct = sum(r["acc"] * r["n_samples"] for r in all_results)
        total_samples = sum(r["n_samples"] for r in all_results)
        weighted_avg = total_correct / total_samples if total_samples > 0 else 0
        logger.info(f"TabPFN: {weighted_avg:.4f}")

   

if __name__ == "__main__":
    main()