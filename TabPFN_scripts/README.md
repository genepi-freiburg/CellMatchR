# CellMatchR – TabPFN

Minimal setup to run TabPFN trained on four reference transcriptomic kidney profile datasets. Reference datasets are hosted as a HuggingFace [dataset](https://huggingface.co/datasets/samuelboehm/cellmatchr) repository.

## What this code does:

1. Downloads reference and test datasets together with a markergene selection automatically from HuggingFace on first run.
2. Selects the intersection of marker genes shared between reference and test data. Marker genes need to be shared between at least one reference and the test dataset.
3. Log-transforms expression values and trains TabPFN.
4. Reports per-dataset and weighted-average accuracy (when true labels are available).

## Setup

```bash
pip install -r requirements.txt
```

## Usage

**Demo mode** — runs on the built-in test datasets:
```bash
python main.py
```

**Your own data**  provide a CSV with cells as rows and genes as columns:
```bash
python main.py --csv /path/to/your/data.csv
```

Gene column headers must be gene symbols (case-insensitive, automatically uppercased). The CSV should not contain any other columns, only raw expression counts in counts per million. When a CSV is provided, the built-in test datasets are skipped and only your data is classified.

Optionally restrict which reference datasets to train on:
```bash
python main.py --csv /path/to/your/data.csv --reference_datasets KPMP Park
```

## Project structure

| File | Description |
|---|---|
| `main.py` | Training, prediction, and evaluation loop |
| `utils.py` | Data loading, gene filtering, HuggingFace download |
| `config.py` | Dataset registry and paths |