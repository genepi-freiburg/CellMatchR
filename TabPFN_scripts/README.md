# CellMatchR – TabPFN

Minimal setup to run TabPFN trained on four reference transcriptomic kidney profile datasets. Reference datasets are hosted as a HuggingFace **dataset** repository.

## What this code does

1. Downloads reference and test datasets together with a marker gene selection automatically from HuggingFace on first run.
2. Selects the intersection of marker genes shared between reference and test data. Marker genes need to be shared between at least one reference and the test dataset.
3. Computes log(CPM) values from raw counts and trains TabPFN.
4. Reports cell type similarity probabilities (figures and tables as csv).
5. Reports per-dataset and weighted-average accuracy (when true labels are available).

## Project structure

| File | Description |
|------|-------------|
| `main.py` | Training, prediction, and evaluation loop |
| `utils.py` | Data loading, gene filtering, HuggingFace download |
| `config.py` | Dataset registry and paths |

---

## Setup

### Windows

<details>
<summary>Click to expand Windows setup instructions</summary>

#### 1. Download scripts
Download the repository scripts to a dedicated folder on your machine.

#### 2. Open PowerShell and start WSL
```powershell
wsl
```

If WSL is not available yet, install it first (requires admin privileges):
```powershell
wsl --install
```
This installs Ubuntu by default (tested with Ubuntu 24.04.4 LTS). A restart may be required.

#### 3. Install pip (if not available)
```bash
sudo apt update && sudo apt install python3-pip -y
```

#### 4. Install venv
```bash
sudo apt install python3.12-venv
```

#### 5. Create a virtual environment
> **Note:** The virtual environment must be created on the Linux filesystem, not on a mounted Windows drive (e.g. `/mnt/d/`), as Windows drives do not support the required Linux file permissions.

```bash
python3 -m venv ~/venvs/cellmatchr
```

#### 6. Activate the virtual environment
```bash
source ~/venvs/cellmatchr/bin/activate
```
You will see `(cellmatchr)` appear in your prompt. **Steps 1–5 only need to be done once.** From now on, just start WSL and run this activation command before using the tool.

#### 7. Install dependencies
Navigate to your scripts folder (Windows drives are accessible under `/mnt/`):
```bash
pip install -r /mnt/c/path/to/your/scripts/requirements.txt
```

#### 8. Set up HuggingFace access token for TabPFN
TabPFN requires a HuggingFace account and access token:
1. Follow the instructions at: https://docs.priorlabs.ai/how-to-access-gated-models
2. Create an access token on HuggingFace
3. Log in via the terminal:
```bash
hf auth login
```
4. Enter your token when prompted

#### 9. Run CellMatchR
```bash
python /mnt/c/path/to/your/scripts/main.py --csv /mnt/c/path/to/your/data.csv
```

</details>

---

### Linux

<details>
<summary>Click to expand Linux setup instructions</summary>

#### 1. Download scripts
Download the repository scripts to a dedicated folder.

#### 2. Create a virtual environment
```bash
python3 -m venv ~/venvs/cellmatchr
```

#### 3. Activate the virtual environment
```bash
source ~/venvs/cellmatchr/bin/activate
```
You will see `(cellmatchr)` appear in your prompt. Run this activation command each time before using the tool.

#### 4. Install dependencies
```bash
pip install -r /path/to/your/scripts/requirements.txt
```

#### 5. Set up HuggingFace access token for TabPFN
TabPFN requires a HuggingFace account and access token:
1. Follow the instructions at: https://docs.priorlabs.ai/how-to-access-gated-models
2. Create an access token on HuggingFace
3. Log in via the terminal:
```bash
hf auth login
```
4. Enter your token when prompted

#### 6. Run CellMatchR
```bash
python /path/to/your/scripts/main.py --csv /path/to/your/data.csv
```

</details>

---

### macOS

<details>
<summary>Click to expand macOS setup instructions</summary>

#### 1. Download scripts
Download the repository scripts to a dedicated folder.

#### 2. Ensure Python 3 is installed
macOS does not ship with Python 3 by default on older versions. Install it via [Homebrew](https://brew.sh) if needed:
```bash
brew install python
```

#### 3. Create a virtual environment
```bash
python3 -m venv ~/venvs/cellmatchr
```

#### 4. Activate the virtual environment
```bash
source ~/venvs/cellmatchr/bin/activate
```
You will see `(cellmatchr)` appear in your prompt. Run this activation command each time before using the tool.

#### 5. Install dependencies
```bash
pip install -r /path/to/your/scripts/requirements.txt
```

#### 6. Set up HuggingFace access token for TabPFN
TabPFN requires a HuggingFace account and access token:
1. Follow the instructions at: https://docs.priorlabs.ai/how-to-access-gated-models
2. Create an access token on HuggingFace
3. Log in via the terminal:
```bash
hf auth login
```
4. Enter your token when prompted

#### 7. Run CellMatchR
```bash
python /path/to/your/scripts/main.py --csv /path/to/your/data.csv
```

</details>

---

## Usage

### Demo mode
Runs on the built-in test datasets:
```bash
python main.py
```

### Your own data
Provide a CSV with cells as rows and genes as columns:
```bash
python main.py --csv /path/to/your/data.csv
```

### Restrict reference datasets
Optionally restrict which reference datasets to train on:
```bash
python main.py --csv /path/to/your/data.csv --reference_datasets KPMP Park
```

---

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

Optionally, add a `meta_target` column with a cell type label per cell to your CSV. If present, an accuracy is calculated and the labels are used as titles in the output probability plots. This is useful for verifying predictions against expected annotations.

> **Note:** Our model is trained on the following cell types: `CD`, `CNT`, `DCT`, `EC`, `ENDO`, `FIB`, `IMM`, `LOH`, `POD`, `PT`. If your `meta_target` column contains other cell type labels, the reported accuracy will not be meaningful, but the labels will still be used for visualization purposes.
