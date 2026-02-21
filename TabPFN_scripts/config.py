import os
from pathlib import Path


REQUIRED_FILES = {
    "KPMP.parquet":             ("reference", None),
    "Park.parquet":             ("reference", None),
    "Ransick.parquet":          ("reference", None),
    "Zhang.parquet":            ("reference", None),
    "Chen.parquet":             ("test", "meta_celltype_tubular"),
    "cortex_medulla.parquet": ("test", "meta_celltype_coarse"),
    "markergenes_all.csv":      ("genelist", None),
}

DATAFOLDER = os.path.join(Path(__file__).parent, 'datasets')

DATA_REPO_ID = "samuelboehm/cellmatchr"
