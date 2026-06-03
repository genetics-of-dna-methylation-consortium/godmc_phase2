# resources/genetics/ld_aggregate.py
"""Pure-Python central LD aggregation for section 15b.

This module deliberately does NOT import hail. Every code path that touches
Hail BlockMatrix IO lives in ld_hail.py; this module operates on NumPy arrays,
pandas DataFrames, and parquet files so its unit tests run without a JVM.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

PRECURSOR_SCHEMA_VERSION = "15b-precursor-v1"
DEFAULT_PANEL_SPEC_VERSION = "0.1.0"
DEFAULT_MAF_THRESHOLD = 0.01
DEFAULT_MIN_ADJ_DIAG = 0.0

VARIANT_COLUMNS = [
    "stable_id", "variant_id", "chr", "pos", "ref", "alt",
    "membership_count", "b_intercept", "b_age", "b_sex",
    "a_diag", "n_nonmissing", "n_imputed",
]

_VARIANT_DTYPES = {
    "stable_id": "int64", "variant_id": "object", "chr": "object",
    "pos": "int64", "ref": "object", "alt": "object",
    "membership_count": "int64", "b_intercept": "float64",
    "b_age": "float64", "b_sex": "float64", "a_diag": "float64",
    "n_nonmissing": "int64", "n_imputed": "int64",
}

PAIR_DTYPE = np.dtype([
    ("pos_i", "<i8"), ("sid_i", "<i8"),
    ("pos_j", "<i8"), ("sid_j", "<i8"),
    ("value", "<f8"),
])
PAIR_KEY_FIELDS = ["pos_i", "sid_i", "pos_j", "sid_j"]


def precursor_paths(precursor_dir: str | Path) -> dict[str, Path]:
    """Return the canonical file/dir paths inside a precursor directory."""
    d = Path(precursor_dir)
    return {
        "root": d,
        "manifest": d / "precursor_manifest.json",
        "d": d / "D.npy",
        "variants": d / "variants.parquet",
        "pairs_dir": d / "A_pairs",
        "lock": d / ".lock",
    }


def _atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    tmp.write_text(text, encoding="utf-8")
    os.replace(tmp, path)


def write_precursor_manifest(precursor_dir: str | Path, manifest: dict) -> None:
    """Atomically write the precursor manifest (temp-then-rename)."""
    path = precursor_paths(precursor_dir)["manifest"]
    _atomic_write_text(path, json.dumps(manifest, indent=2, sort_keys=True) + "\n")


def read_precursor_manifest(precursor_dir: str | Path) -> dict | None:
    """Read the precursor manifest, or None if it does not exist yet."""
    path = precursor_paths(precursor_dir)["manifest"]
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _empty_variant_table() -> pd.DataFrame:
    return pd.DataFrame({c: pd.Series(dtype=_VARIANT_DTYPES[c]) for c in VARIANT_COLUMNS})


def read_variant_table(precursor_dir: str | Path) -> pd.DataFrame:
    """Read the append-only master variant table, or an empty typed frame."""
    path = precursor_paths(precursor_dir)["variants"]
    if not path.is_file():
        return _empty_variant_table()
    return pd.read_parquet(path)[VARIANT_COLUMNS].astype(_VARIANT_DTYPES)


def write_variant_table(precursor_dir: str | Path, df: pd.DataFrame) -> None:
    """Atomically write the master variant table as parquet."""
    path = precursor_paths(precursor_dir)["variants"]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")
    df[VARIANT_COLUMNS].astype(_VARIANT_DTYPES).to_parquet(tmp, index=False)
    os.replace(tmp, path)
