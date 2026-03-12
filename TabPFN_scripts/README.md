# CellMatchR – TabPFN

Minimal setup to run TabPFN trained on four reference transcriptomic kidney profile datasets. Reference datasets are hosted as a HuggingFace dataset repository.

## What this code does

1. Downloads reference and test datasets together with a marker gene selection automatically from HuggingFace on first run.
2. Selects the intersection of marker genes shared between reference and test data. Marker genes need to be shared between at least one reference and the test dataset.
3. Computes log(CPM + 1) values from raw counts and trains TabPFN.
4. Reports cell type similarity probabilities (figures and tables as csv).
5. Reports per-dataset and weighted-average accuracy (when true labels are available).

## Project structure

| File | Description |
|------|-------------|
| `main.py` | Training, prediction, and evaluation loop |
| `utils.py` | Data loading, gene filtering, HuggingFace download |
| `config.py` | Dataset registry and paths |

---

## Quick Start with Google Colab

For the easiest experience, use our Google Colab notebook:

🔗 **[Open in Colab](https://colab.research.google.com/github/genepi-freiburg/CellMatchR/blob/main/TabPFN_scripts/CellMatchR_TabPFN.ipynb)**

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/genepi-freiburg/CellMatchR/blob/main/TabPFN_scripts/CellMatchR_TabPFN.ipynb)

The Colab notebook provides:
- Pre-configured environment with all dependencies
- Free GPU acceleration
- Step-by-step guide through the CellMatchR workflow
- No local installation required

---

## Local Setup

**Prerequisites:** Python 3.8+, pip, basic command line knowledge

### Quick Setup

```bash
# Clone repository
git clone https://github.com/genepi-freiburg/CellMatchR.git
cd CellMatchR/TabPFN_scripts

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # Linux/macOS
# venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt

# Authenticate with HuggingFace (required for TabPFN)
huggingface-cli login

# Run CellMatchR
python main.py --csv path/to/your/data.csv
```

### HuggingFace Authentication

TabPFN requires a HuggingFace account to downlad the pre-trained model:
1. Get token: https://huggingface.co/settings/tokens
2. Authenticate: `huggingface-cli login`
3. Paste your token when prompted


## Usage Examples

```bash
# Demo with built-in test data
python main.py

# Your own CSV data
python main.py --csv path/to/data.csv

# Specific reference datasets
python main.py --csv data.csv --reference_datasets KPMP Park

# Help
python main.py --help
```

---

### Runtime
Runtime varies depending on hardware and dataset size. If a GPU is available, predictions are typically fast (seconds to low minutes). When running on CPU only, runtimes of >10 minutes are possible for larger datasets. The tool automatically uses a GPU if detected, and falls back to CPU otherwise.

## Input data format

- One **row per cell**, one **column per gene**
- Gene column headers must be **gene symbols** (case-insensitive, automatically uppercased)
- The CSV should contain only **raw expression counts**
- When a CSV is provided, the built-in test datasets are skipped and only your data is classified

### Example input (without labels)

| Slc12a1 | Umod | Nphs1 | Cdh1 | Lrp2 |
|---------|------|-------|------|------|
| 0 | 0 | 0 | 524 | 1203 |
| 0 | 312 | 0 | 0 | 0 |
| 87 | 0 | 0 | 0 | 0 |
| 0 | 0 | 445 | 0 | 0 |

### Example input (with optional `meta_target` label column)

| meta_target | Slc12a1 | Umod | Nphs1 | Cdh1 | Lrp2 |
|-------------|---------|------|-------|------|------|
| PT | 0 | 0 | 0 | 524 | 1203 |
| LOH | 0 | 312 | 0 | 0 | 0 |
| LOH | 87 | 0 | 0 | 0 | 0 |
| POD | 0 | 0 | 445 | 0 | 0 |

The optional `meta_target` is recomended. If present, the labels are used as titles in the output probability plots. This is useful for verifying predictions against expected annotations. If the given celltypes match the training celltype annotations, an accuracy value of the matching is reported.

> **Note:** Our model is trained on the following cell types: `CD`, `CNT`, `DCT`, `EC`, `ENDO`, `FIB`, `IMM`, `LOH`, `POD`, `PT`. If your `meta_target` column contains other cell type labels, the reported accuracy will not be meaningful.