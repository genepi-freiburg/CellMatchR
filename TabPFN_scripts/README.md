# CellMatchR – TabPFN

Minimal setup to run TabPFN trained on four reference transcriptomic kidney profile datasets. Reference datasets are hosted as a HuggingFace [dataset](https://huggingface.co/datasets/samuelboehm/cellmatchr) repository.

## What this code does:

1. Downloads reference and test datasets together with a markergene selection automatically from HuggingFace on first run.
2. Selects the intersection of marker genes shared between reference and test data. Marker genes need to be shared between at least one reference and the test dataset.
3. Log-transforms expression values and trains TabPFN.
4. Reports per-dataset and weighted-average accuracy.

## Setup

```bash
pip install -r requirements.txt
```

## Usage

```bash
python main.py
```

Datasets are cached locally in `./datasets/` after the first run.

## Project structure

| File | Description |
|---|---|
| `main.py` | Training, prediction, and evaluation loop |
| `utils.py` | Data loading, gene filtering, HuggingFace download |
| `config.py` | Dataset registry and paths |