"""Project paths shared by comparison scripts.

The repository root is detected from this file. Set BGROOT only when running
against a copy located elsewhere.
"""
from pathlib import Path
import os

ROOT = Path(os.environ.get("BGROOT", Path(__file__).resolve().parents[2])).resolve()
COMPARISON = ROOT / "comparison"
TABLES = COMPARISON / "tables"
FIGURES = COMPARISON / "figures"
DATASETS = ROOT / "datasets"
ZHUANG = DATASETS / "zhang_zhuang_2023_merfish"
HAO = DATASETS / "hao_2024_stereo_seq"
SUYAMA = DATASETS / "suyama_2026_xenium"
SILETTI = DATASETS / "siletti_2023_snrna"
KOZAREVA = DATASETS / "kozareva_2021_snrna"
