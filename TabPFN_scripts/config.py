import os
from pathlib import Path


REQUIRED_FILES = {
    "KPMP.parquet" : "reference",
    "Park.parquet" : "reference",
    "Ransick.parquet" : "reference",
    "Zhang.parquet" : "reference",
    "Chen.parquet" : "test",
    "coarse_celltypes.parquet" : "test",
    "markergenes_all.csv" : "genelist"
}

DATAFOLDER = os.path.join(Path(__file__).parent, 'datasets')

DATA_REPO_ID = "samuelboehm/cellmatchr"
