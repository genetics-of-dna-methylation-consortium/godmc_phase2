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


CONTRACT_FIELDS = [
    "genome_build", "schema_id", "matrix_columns", "sex_factor_recode",
    "variants_schema_version", "radius_bp", "block_size",
]


def extract_contract(cohort_manifest: dict) -> dict:
    """Pull the cross-cohort contract fields out of a 15a cohort manifest."""
    cov = cohort_manifest["covariate_schema"]
    a = cohort_manifest["A_blocks"]
    return {
        "genome_build": cohort_manifest["genome_build"],
        "schema_id": cov["schema_id"],
        "matrix_columns": list(cov["matrix_columns"]),
        "sex_factor_recode": dict(cov["sex_factor_recode"]),
        "variants_schema_version": cohort_manifest["variant_index"]["schema_version"],
        "radius_bp": a["radius_bp"],
        "block_size": a["block_size"],
    }


def validate_contract(expected: dict, candidate: dict) -> None:
    """Hard-fail if a cohort's contract disagrees with the precursor's."""
    for field in CONTRACT_FIELDS:
        if expected[field] != candidate[field]:
            raise ValueError(
                f"Cohort contract mismatch on '{field}': precursor has "
                f"{expected[field]!r}, cohort has {candidate[field]!r}. "
                "All cohorts in a pooled panel must share this value."
            )


_COHORT_VARIANT_DTYPES = {
    "chr": "object", "pos": "int64", "ref": "object", "alt": "object",
    "variant_id": "object", "n_nonmissing": "int64", "n_imputed": "int64",
    "genotype_mean": "float64",
}


def read_cohort_assets(cohort_dir: str | Path) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    """Read a 15a cohort's variants.tsv.gz, B.npy, and D.npy.

    Returns (variants_df in file/row order, B (n_variants x 3), D (3 x 3)).
    Hard-fails if B's row count does not match the variant count.
    """
    cohort_dir = Path(cohort_dir)
    variants = pd.read_csv(
        cohort_dir / "variants.tsv.gz", sep="\t", dtype=_COHORT_VARIANT_DTYPES,
    )
    b_mat = np.load(cohort_dir / "B.npy", allow_pickle=False).astype(np.float64)
    d_mat = np.load(cohort_dir / "D.npy", allow_pickle=False).astype(np.float64)
    if b_mat.shape[0] != len(variants):
        raise ValueError(
            f"B.npy has {b_mat.shape[0]} rows but variants.tsv.gz has "
            f"{len(variants)} rows; they must align one-to-one"
        )
    if b_mat.shape[1] != 3 or d_mat.shape != (3, 3):
        raise ValueError(
            f"Expected B (n x 3) and D (3 x 3); got B{b_mat.shape}, D{d_mat.shape}"
        )
    return variants, b_mat, d_mat


def merge_variant_table(
    table: pd.DataFrame,
    cohort_variants: pd.DataFrame,
    b_matrix: np.ndarray,
    next_stable_id: int,
) -> tuple[pd.DataFrame, dict[str, int], int]:
    """Add one cohort's per-variant data into the master table.

    Assigns a fresh stable_id to each unseen variant_id, sums the three B
    columns, n_nonmissing, n_imputed, and bumps membership_count. a_diag is
    untouched here (filled by the A-merge). Returns (updated table,
    {variant_id: stable_id} for this cohort, next free stable_id).
    """
    existing = table.set_index("variant_id")
    known_ids = dict(zip(existing.index, existing["stable_id"]))

    id_map: dict[str, int] = {}
    new_rows: list[dict] = []
    for pos_in_cohort, row in enumerate(cohort_variants.itertuples(index=False)):
        vid = row.variant_id
        b = b_matrix[pos_in_cohort]
        if vid in known_ids:
            sid = int(known_ids[vid])
            existing.loc[vid, "b_intercept"] += b[0]
            existing.loc[vid, "b_age"] += b[1]
            existing.loc[vid, "b_sex"] += b[2]
            existing.loc[vid, "n_nonmissing"] += int(row.n_nonmissing)
            existing.loc[vid, "n_imputed"] += int(row.n_imputed)
            existing.loc[vid, "membership_count"] += 1
        else:
            sid = next_stable_id
            next_stable_id += 1
            known_ids[vid] = sid
            new_rows.append({
                "stable_id": sid, "variant_id": vid, "chr": row.chr,
                "pos": int(row.pos), "ref": row.ref, "alt": row.alt,
                "membership_count": 1, "b_intercept": float(b[0]),
                "b_age": float(b[1]), "b_sex": float(b[2]), "a_diag": 0.0,
                "n_nonmissing": int(row.n_nonmissing), "n_imputed": int(row.n_imputed),
            })
        id_map[vid] = sid

    updated = existing.reset_index()[VARIANT_COLUMNS]
    if new_rows:
        updated = pd.concat([updated, pd.DataFrame(new_rows)[VARIANT_COLUMNS]],
                            ignore_index=True)
    return updated.astype(_VARIANT_DTYPES), id_map, next_stable_id


def window_stops(positions: np.ndarray, radius_bp: int) -> np.ndarray:
    """Exclusive column stop per row: first index with pos > pos_i + radius."""
    return np.searchsorted(positions, positions + radius_bp, side="right")


def extract_chunk_entries(
    dense: np.ndarray,
    row_start: int,
    col_start: int,
    positions: np.ndarray,
    sids: np.ndarray,
    radius_bp: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split a dense A chunk into diagonal updates and off-diagonal pair rows.

    `positions`/`sids` are whole-chromosome arrays (canonical order). `dense`
    covers chromosome rows [row_start, row_start+n_rows) and columns starting
    at chromosome index `col_start`. Returns (diag_sids, diag_vals, offdiag)
    where offdiag is a structured array of PAIR_DTYPE with pos_i <= pos_j.
    """
    stops = window_stops(positions, radius_bp)
    n_rows = dense.shape[0]
    diag_sids = np.empty(n_rows, dtype=np.int64)
    diag_vals = np.empty(n_rows, dtype=np.float64)
    off_chunks: list[np.ndarray] = []

    for r in range(n_rows):
        gi = row_start + r
        diag_sids[r] = sids[gi]
        diag_vals[r] = dense[r, gi - col_start]
        j0, j1 = gi + 1, int(stops[gi])
        if j1 <= j0:
            continue
        cols = np.arange(j0, j1)
        rec = np.empty(j1 - j0, dtype=PAIR_DTYPE)
        rec["pos_i"] = positions[gi]
        rec["sid_i"] = sids[gi]
        rec["pos_j"] = positions[cols]
        rec["sid_j"] = sids[cols]
        rec["value"] = dense[r, cols - col_start]
        off_chunks.append(rec)

    offdiag = (np.concatenate(off_chunks) if off_chunks
               else np.empty(0, dtype=PAIR_DTYPE))
    return diag_sids, diag_vals, offdiag
