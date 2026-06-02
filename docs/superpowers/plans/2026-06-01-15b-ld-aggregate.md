# Section 15b Central LD Aggregation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the pre-Hail 15b scaffold with a two-phase central aggregator — `accumulate` (add one cohort's `A`/`B`/`D` into a mutable additive precursor) and `finalise` (resolve the cohort intersection, covariate-adjust, and emit an immutable pooled LD correlation panel `R`).

**Architecture:** Per-variant additive statistics (`B` row, `A` diagonal, frequency inputs, cohort-membership count) live in an append-only `variants.parquet`; the only large structure is an off-diagonal windowed `A` pair store, sorted by genomic position so `accumulate` is a linear sorted-merge and `finalise` streams it in pooled-row order. All aggregation and the rank-3 covariate adjustment run as bounded per-chromosome NumPy chunks (never a dense p×p); Hail is used only to read back cohort `A_blocks` and write float32 `R_blocks`.

**Tech Stack:** Python 3.12, NumPy 2.3, SciPy 1.16, pandas 2.x, pyarrow 24 (parquet), Hail 0.2.135 (BlockMatrix IO only), pytest. Spec: `docs/superpowers/specs/2026-06-01-15b-ld-aggregate-design.md`.

**Conventions used throughout:**
- Covariate matrix columns are `[intercept, Age_numeric, Sex_factor]`, so `B[:,0]=b_intercept`, `B[:,1]=b_age`, `B[:,2]=b_sex`.
- Tests live in the git-excluded `tests/genetics/` tree (same `hail_env`, same `conftest.py` that injects `resources/genetics` onto `sys.path`). Hail tests use `pytest.importorskip("hail")` and/or the `hail_session` fixture.
- Run tests with: `conda run -n hail_env python -m pytest tests/genetics/<file>::<test> -v` (substitute the project's usual pytest invocation if different).
- Pure-module tests (no Hail) can run under any env with numpy/pandas/pyarrow.

---

## Shared definitions (referenced by every task)

These names are defined in Task 1 and reused verbatim later. They are listed here so tasks read out of order stay consistent.

```python
# resources/genetics/ld_aggregate.py — module constants
PRECURSOR_SCHEMA_VERSION = "15b-precursor-v1"
DEFAULT_PANEL_SPEC_VERSION = "0.1.0"
DEFAULT_MAF_THRESHOLD = 0.01
DEFAULT_MIN_ADJ_DIAG = 0.0  # drop variants whose adjusted diagonal is <= this

VARIANT_COLUMNS = [
    "stable_id", "variant_id", "chr", "pos", "ref", "alt",
    "membership_count", "b_intercept", "b_age", "b_sex",
    "a_diag", "n_nonmissing", "n_imputed",
]

import numpy as np

PAIR_DTYPE = np.dtype([
    ("pos_i", "<i8"), ("sid_i", "<i8"),
    ("pos_j", "<i8"), ("sid_j", "<i8"),
    ("value", "<f8"),
])
PAIR_KEY_FIELDS = ["pos_i", "sid_i", "pos_j", "sid_j"]
```

---

## Task 1: Precursor scaffolding — paths, manifest IO, table schemas

**Files:**
- Create: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/genetics/test_ld_aggregate.py
import json
import numpy as np
import pandas as pd
import pytest

import ld_aggregate as agg


def test_precursor_paths_are_under_dir(tmp_path):
    p = agg.precursor_paths(tmp_path)
    assert p["manifest"] == tmp_path / "precursor_manifest.json"
    assert p["d"] == tmp_path / "D.npy"
    assert p["variants"] == tmp_path / "variants.parquet"
    assert p["pairs_dir"] == tmp_path / "A_pairs"


def test_manifest_round_trip_is_atomic(tmp_path):
    manifest = {"schema_version": agg.PRECURSOR_SCHEMA_VERSION, "n_cohorts": 2}
    agg.write_precursor_manifest(tmp_path, manifest)
    assert agg.read_precursor_manifest(tmp_path) == manifest
    # no leftover temp files
    assert not list(tmp_path.glob("*.tmp*"))


def test_read_manifest_missing_returns_none(tmp_path):
    assert agg.read_precursor_manifest(tmp_path) is None


def test_empty_variant_table_has_schema(tmp_path):
    df = agg.read_variant_table(tmp_path)
    assert list(df.columns) == agg.VARIANT_COLUMNS
    assert len(df) == 0


def test_variant_table_round_trip(tmp_path):
    df = pd.DataFrame(
        [[0, "1:100:G:A", "1", 100, "G", "A", 1, 3.0, 90.0, 1.5, 5.0, 4, 0]],
        columns=agg.VARIANT_COLUMNS,
    )
    agg.write_variant_table(tmp_path, df)
    back = agg.read_variant_table(tmp_path)
    pd.testing.assert_frame_equal(back, df)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'ld_aggregate'`.

- [ ] **Step 3: Write minimal implementation**

```python
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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -v`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py
git commit -m "feat(15b): precursor paths, manifest IO, and table schemas"
```

---

## Task 2: Contract extraction and validation

**Files:**
- Modify: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

- [ ] **Step 1: Write the failing test**

```python
def _cohort_manifest():
    return {
        "study_name": "cohortA",
        "genome_build": "GRCh37",
        "covariate_schema": {
            "schema_id": "intercept_age_sex",
            "matrix_columns": ["intercept", "Age_numeric", "Sex_factor"],
            "sex_factor_recode": {"M": 1.0, "F": 2.0},
        },
        "variant_index": {"schema_version": "v0.3-with-genotype-stats"},
        "A_blocks": {"radius_bp": 1_000_000, "block_size": 4096},
    }


def test_extract_contract_pulls_expected_fields():
    c = agg.extract_contract(_cohort_manifest())
    assert c == {
        "genome_build": "GRCh37",
        "schema_id": "intercept_age_sex",
        "matrix_columns": ["intercept", "Age_numeric", "Sex_factor"],
        "sex_factor_recode": {"M": 1.0, "F": 2.0},
        "variants_schema_version": "v0.3-with-genotype-stats",
        "radius_bp": 1_000_000,
        "block_size": 4096,
    }


def test_validate_contract_matches():
    c = agg.extract_contract(_cohort_manifest())
    agg.validate_contract(c, c)  # no raise


def test_validate_contract_mismatch_raises():
    base = agg.extract_contract(_cohort_manifest())
    other = dict(base, block_size=1024)
    with pytest.raises(ValueError, match="block_size"):
        agg.validate_contract(base, other)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k contract -v`
Expected: FAIL with `AttributeError: module 'ld_aggregate' has no attribute 'extract_contract'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k contract -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): cohort contract extraction and validation"
```

---

## Task 3: Read a cohort's variants + B + D

**Files:**
- Modify: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

The 15a `variants.tsv.gz` columns are `chr,pos,ref,alt,variant_id,n_nonmissing,n_imputed,genotype_mean`; `B.npy` rows align to that order; `D.npy` is 3×3.

- [ ] **Step 1: Write the failing test**

```python
import gzip


def _write_cohort_dir(tmp_path, variant_rows, b, d):
    out = tmp_path / "cohort"
    out.mkdir()
    header = "chr\tpos\tref\talt\tvariant_id\tn_nonmissing\tn_imputed\tgenotype_mean\n"
    body = "".join("\t".join(map(str, r)) + "\n" for r in variant_rows)
    with gzip.open(out / "variants.tsv.gz", "wt") as fh:
        fh.write(header + body)
    np.save(out / "B.npy", np.asarray(b, dtype=np.float64))
    np.save(out / "D.npy", np.asarray(d, dtype=np.float64))
    return out


def test_read_cohort_assets(tmp_path):
    rows = [
        ("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0),
        ("1", 200, "C", "T", "1:200:C:T", 3, 1, 0.5),
    ]
    b = [[4.0, 90.0, 1.5], [2.0, 70.0, 1.0]]
    d = [[4, 180, 6], [180, 8200, 270], [6, 270, 10]]
    out = _write_cohort_dir(tmp_path, rows, b, d)

    variants, b_mat, d_mat = agg.read_cohort_assets(out)
    assert list(variants["variant_id"]) == ["1:100:G:A", "1:200:C:T"]
    assert variants["pos"].tolist() == [100, 200]
    np.testing.assert_array_equal(b_mat, np.asarray(b))
    np.testing.assert_array_equal(d_mat, np.asarray(d))
    assert b_mat.shape == (2, 3)


def test_read_cohort_assets_rejects_b_row_mismatch(tmp_path):
    rows = [("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0)]
    out = _write_cohort_dir(tmp_path, rows, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
                            np.eye(3))
    with pytest.raises(ValueError, match="rows"):
        agg.read_cohort_assets(out)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k cohort_assets -v`
Expected: FAIL with `AttributeError: ... 'read_cohort_assets'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k cohort_assets -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): read cohort variants/B/D assets"
```

---

## Task 4: Per-variant merge into the master table

**Files:**
- Modify: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

Merge adds the three `B` columns, bumps `membership_count`, adds `n_nonmissing`/`n_imputed`, and assigns a new `stable_id` to any unseen `variant_id`. `a_diag` is initialised to 0 here and filled by the A-merge (Task 8). Returns the updated table and a `variant_id -> stable_id` map for the cohort's variants.

- [ ] **Step 1: Write the failing test**

```python
def _variants_df(rows):
    cols = ["chr", "pos", "ref", "alt", "variant_id", "n_nonmissing",
            "n_imputed", "genotype_mean"]
    return pd.DataFrame(rows, columns=cols)


def test_merge_variant_table_assigns_ids_and_sums():
    table = agg._empty_variant_table()
    cohort = _variants_df([
        ("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0),
        ("1", 200, "C", "T", "1:200:C:T", 3, 1, 0.5),
    ])
    b = np.array([[4.0, 90.0, 1.5], [2.0, 70.0, 1.0]])

    table, id_map, next_id = agg.merge_variant_table(table, cohort, b, next_stable_id=0)
    assert next_id == 2
    assert id_map == {"1:100:G:A": 0, "1:200:C:T": 1}
    assert table.loc[table.variant_id == "1:100:G:A", "b_intercept"].item() == 4.0
    assert table["membership_count"].tolist() == [1, 1]

    # second cohort: one shared variant, one new
    cohort2 = _variants_df([
        ("1", 100, "G", "A", "1:100:G:A", 5, 0, 1.2),
        ("1", 300, "A", "G", "1:300:A:G", 5, 0, 0.8),
    ])
    b2 = np.array([[6.0, 100.0, 2.0], [3.0, 80.0, 1.0]])
    table, id_map2, next_id = agg.merge_variant_table(table, cohort2, b2, next_stable_id=next_id)

    assert next_id == 3
    assert id_map2 == {"1:100:G:A": 0, "1:300:A:G": 2}
    shared = table.set_index("variant_id").loc["1:100:G:A"]
    assert shared["membership_count"] == 2
    assert shared["b_intercept"] == 10.0       # 4 + 6
    assert shared["n_nonmissing"] == 9         # 4 + 5
    assert table.set_index("variant_id").loc["1:300:A:G"]["membership_count"] == 1
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k merge_variant_table -v`
Expected: FAIL with `AttributeError: ... 'merge_variant_table'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k merge_variant_table -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): per-variant merge into master table"
```

---

## Task 5: Extract diagonal + off-diagonal entries from one A chunk

**Files:**
- Modify: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

A 15a chunk is a dense `(n_rows, n_cols)` block where row `r` (chromosome index `row_start+r`) holds the 1 Mb window `cols [global_i, stop)` in chromosome-local-to-chunk coordinates (`col_start == row_start`). This is pure NumPy: given the dense block, positions, and per-row stable_ids, split into diagonal updates and off-diagonal `(pos_i, sid_i, pos_j, sid_j, value)` rows. Window membership is recomputed from positions (robust to legitimate zero values), matching 15a's `stops`.

- [ ] **Step 1: Write the failing test**

```python
def test_extract_chunk_entries_splits_diag_and_offdiag():
    # 3 variants on one chrom, positions 100/200/2_000_000, radius 1Mb.
    # pos 100 windows {100,200}; pos 200 windows {200,2.0e6}? 200+1e6 < 2e6 so no.
    positions = np.array([100, 200, 2_000_000], dtype=np.int64)
    sids = np.array([10, 11, 12], dtype=np.int64)
    radius = 1_000_000
    # full chunk: rows 0..3, cols 0..3 (col_start=0)
    dense = np.array([
        [5.0, 2.0, 0.0],   # row0: diag=5 (v0,v0), offdiag (v0,v1)=2
        [0.0, 6.0, 0.0],   # row1: diag=6 (v1,v1); (v1,v2) out of window
        [0.0, 0.0, 7.0],   # row2: diag=7 (v2,v2)
    ])
    diag_sid, diag_val, offdiag = agg.extract_chunk_entries(
        dense, row_start=0, col_start=0, positions=positions, sids=sids,
        radius_bp=radius,
    )
    assert diag_sid.tolist() == [10, 11, 12]
    assert diag_val.tolist() == [5.0, 6.0, 7.0]
    assert offdiag.shape[0] == 1
    rec = offdiag[0]
    assert (rec["pos_i"], rec["sid_i"], rec["pos_j"], rec["sid_j"], rec["value"]) == (
        100, 10, 200, 11, 2.0)


def test_extract_chunk_entries_offset_chunk():
    # rows 1..3 only, col_start=1 (a later chunk)
    positions = np.array([100, 200, 300], dtype=np.int64)
    sids = np.array([10, 11, 12], dtype=np.int64)
    dense = np.array([
        [6.0, 3.0],   # row global1: diag(v1)=6, (v1,v2)=3
        [0.0, 7.0],   # row global2: diag(v2)=7
    ])
    diag_sid, diag_val, offdiag = agg.extract_chunk_entries(
        dense, row_start=1, col_start=1, positions=positions, sids=sids,
        radius_bp=1_000_000,
    )
    assert diag_sid.tolist() == [11, 12]
    assert offdiag[0]["sid_i"] == 11 and offdiag[0]["sid_j"] == 12
    assert offdiag[0]["value"] == 3.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k extract_chunk -v`
Expected: FAIL with `AttributeError: ... 'extract_chunk_entries'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py

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
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k extract_chunk -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): extract diagonal/off-diagonal entries from an A chunk"
```

---

## Task 6: Streaming sorted-merge-add of the pair store

**Files:**
- Modify: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

`sorted_merge_add` merges two key-sorted `PAIR_DTYPE` arrays, summing `value` on equal keys. `merge_pairs_into_file` streams the existing per-chromosome parquet in batches and merges an in-memory incoming sorted array, writing a new sorted parquet atomically.

> **Release-1 scope note (not a placeholder):** the incoming cohort-chromosome off-diagonal array is held in memory while the *existing precursor* file (the side that grows with cohort count) is streamed in batches. This is sufficient for the two-cohort pilot. Batching the incoming side as well is a tractability follow-up tied to open question #92 (genome-wide vs fine-mapping scope) and is explicitly out of release-1 scope.

- [ ] **Step 1: Write the failing test**

```python
import pyarrow.parquet as pq


def _pairs(rows):
    arr = np.empty(len(rows), dtype=agg.PAIR_DTYPE)
    for i, r in enumerate(rows):
        arr[i] = r
    return arr


def test_sorted_merge_add_sums_equal_keys():
    a = _pairs([(100, 0, 200, 1, 2.0), (100, 0, 300, 2, 5.0)])
    b = _pairs([(100, 0, 200, 1, 3.0), (150, 3, 400, 4, 1.0)])
    out = agg.sorted_merge_add(a, b)
    assert out["value"].tolist() == [5.0, 5.0, 1.0]      # (100,0,200,1) summed
    assert out["pos_i"].tolist() == [100, 100, 150]
    assert out["pos_j"].tolist() == [200, 300, 400]


def test_merge_pairs_into_file_round_trip(tmp_path):
    path = tmp_path / "chr1.parquet"
    first = _pairs([(100, 0, 200, 1, 2.0)])
    agg.merge_pairs_into_file(path, first, batch_rows=8)
    second = _pairs([(100, 0, 200, 1, 3.0), (300, 2, 400, 3, 1.0)])
    agg.merge_pairs_into_file(path, second, batch_rows=8)

    back = pq.read_table(path).to_pandas()
    assert back["value"].tolist() == [5.0, 1.0]
    assert back["pos_i"].tolist() == [100, 300]
    assert not list(tmp_path.glob("*.tmp*"))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k "merge_add or merge_pairs" -v`
Expected: FAIL with `AttributeError: ... 'sorted_merge_add'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py
import pyarrow as pa
import pyarrow.parquet as pq

_PAIR_ARROW_SCHEMA = pa.schema([
    ("pos_i", pa.int64()), ("sid_i", pa.int64()),
    ("pos_j", pa.int64()), ("sid_j", pa.int64()), ("value", pa.float64()),
])


def _pair_key(arr: np.ndarray) -> np.ndarray:
    """Structured view of the 4 key fields for lexicographic comparison."""
    return arr[PAIR_KEY_FIELDS]


def sorted_merge_add(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Merge two key-sorted PAIR_DTYPE arrays, summing value on equal keys."""
    merged = np.concatenate([a, b])
    order = np.argsort(_pair_key(merged), kind="stable")
    merged = merged[order]
    if merged.size == 0:
        return merged
    keys = _pair_key(merged)
    # boundary[k] True when row k starts a new key group
    same = (keys[1:] == keys[:-1])
    group = np.empty(merged.size, dtype=np.int64)
    group[0] = 0
    np.cumsum(~same, out=group[1:])
    n_groups = int(group[-1]) + 1
    out = merged[np.searchsorted(group, np.arange(n_groups))]
    summed = np.zeros(n_groups, dtype=np.float64)
    np.add.at(summed, group, merged["value"])
    out["value"] = summed
    return out


def _arrow_from_pairs(arr: np.ndarray) -> "pa.Table":
    return pa.table({f: arr[f] for f in _PAIR_ARROW_SCHEMA.names},
                    schema=_PAIR_ARROW_SCHEMA)


def _pairs_from_arrow(batch) -> np.ndarray:
    cols = {name: batch.column(name).to_numpy(zero_copy_only=False)
            for name in _PAIR_ARROW_SCHEMA.names}
    arr = np.empty(batch.num_rows, dtype=PAIR_DTYPE)
    for name in _PAIR_ARROW_SCHEMA.names:
        arr[name] = cols[name]
    return arr


def merge_pairs_into_file(path: str | Path, incoming: np.ndarray, batch_rows: int) -> None:
    """Sorted-merge-add an in-memory incoming pair array into a parquet file.

    Streams the existing file in row batches (the side that grows with cohort
    count); the incoming array is assumed key-sorted and held in memory.
    Writes atomically (temp-then-rename).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    incoming = incoming[np.argsort(_pair_key(incoming), kind="stable")]
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}")

    writer = pq.ParquetWriter(tmp, _PAIR_ARROW_SCHEMA)
    try:
        if not path.is_file():
            writer.write_table(_arrow_from_pairs(sorted_merge_add(
                np.empty(0, dtype=PAIR_DTYPE), incoming)))
        else:
            pf = pq.ParquetFile(path)
            in_pos = 0
            for batch in pf.iter_batches(batch_size=batch_rows):
                existing = _pairs_from_arrow(batch)
                # take incoming rows whose key < last existing key in this batch
                last_key = existing[-1:][PAIR_KEY_FIELDS]
                take = in_pos
                inc_keys = incoming[PAIR_KEY_FIELDS]
                while take < incoming.size and tuple(inc_keys[take]) <= tuple(last_key[0]):
                    take += 1
                chunk_inc = incoming[in_pos:take]
                in_pos = take
                writer.write_table(_arrow_from_pairs(
                    sorted_merge_add(existing, chunk_inc)))
            if in_pos < incoming.size:
                writer.write_table(_arrow_from_pairs(incoming[in_pos:]))
    finally:
        writer.close()
    os.replace(tmp, path)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k "merge_add or merge_pairs" -v`
Expected: PASS.

> Note for the implementer: the batch-merge above sums equal keys *within* the spliced (existing-batch + incoming-slice) window. Because both sides are globally key-sorted and the incoming slice is bounded by the batch's last key, a key shared across the existing/incoming boundary always lands in the same `sorted_merge_add` call. Verify this with the round-trip test before moving on.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): streaming sorted-merge-add for the pair store"
```

---

## Task 7: Hail BlockMatrix read/write adapters

**Files:**
- Modify: `resources/genetics/ld_hail.py`
- Test: `tests/genetics/test_ld_hail.py`

Two thin functions: read a written `A_blocks` chunk back to a dense NumPy array, and write a dense `R` block (sparsified to upper-triangular row intervals) as a float32 BlockMatrix directory mirroring `compute_a_block_banded`'s output.

- [ ] **Step 1: Write the failing test**

```python
# Append to tests/genetics/test_ld_hail.py
import numpy as np
from hail.linalg import BlockMatrix
import ld_hail


def test_read_a_block_chunk_round_trip(hail_session, tmp_path):
    arr = np.array([[5.0, 2.0, 0.0], [0.0, 6.0, 3.0], [0.0, 0.0, 7.0]])
    bm = BlockMatrix.from_numpy(arr, block_size=8)
    out = tmp_path / "chunk"
    bm.write(str(out), overwrite=True)
    back = ld_hail.read_a_block_chunk(out)
    np.testing.assert_array_equal(back, arr)


def test_write_r_block_chunk_sparsifies_row_intervals(hail_session, tmp_path):
    dense = np.array([[1.0, 0.4, 0.2], [0.0, 1.0, 0.5], [0.0, 0.0, 1.0]],
                     dtype=np.float64)
    starts = np.array([0, 1, 2], dtype=np.int64)
    stops = np.array([3, 3, 3], dtype=np.int64)   # full upper triangle within window
    out = tmp_path / "r_chunk"
    ld_hail.write_r_block_chunk(dense, starts, stops, out, block_size=8)
    back = BlockMatrix.read(str(out)).to_numpy()
    assert back.dtype == np.float32
    np.testing.assert_allclose(np.triu(back), np.triu(dense), rtol=1e-6)
    # below-diagonal entries were sparsified away
    assert back[1, 0] == 0.0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_hail.py -k "a_block_chunk or r_block_chunk" -v`
Expected: FAIL with `AttributeError: module 'ld_hail' has no attribute 'read_a_block_chunk'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_hail.py (BlockMatrix already imported at top)

def read_a_block_chunk(chunk_dir: str | Path) -> np.ndarray:
    """Read a written A_blocks chunk BlockMatrix back as a dense float64 array."""
    return BlockMatrix.read(str(chunk_dir)).to_numpy().astype(np.float64, copy=False)


def write_r_block_chunk(
    dense: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
    out_dir: str | Path,
    block_size: int = DEFAULT_A_BLOCK_SIZE,
) -> None:
    """Write a dense R block as a float32 row-interval BlockMatrix directory.

    `starts`/`stops` are per-row chunk-local column intervals (same convention
    as compute_a_block_banded's sparsify_row_intervals call).
    """
    bm = BlockMatrix.from_numpy(dense.astype(np.float32), block_size=block_size)
    bm = bm.sparsify_row_intervals(
        starts=[int(s) for s in starts],
        stops=[int(s) for s in stops],
        blocks_only=False,
    )
    bm.write(str(out_dir), overwrite=True)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_hail.py -k "a_block_chunk or r_block_chunk" -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_hail.py tests/genetics/test_ld_hail.py
git commit -m "feat(15b): Hail A-chunk read and R-chunk write adapters"
```

---

## Task 8: `accumulate` orchestration + CLI + lockfile + atomic commit

**Files:**
- Modify: `resources/genetics/ld_aggregate.py` (orchestration, Hail injected as a callback)
- Rewrite: `resources/genetics/ld_aggregate_stats.py` (CLI `--mode accumulate`)
- Test: `tests/genetics/test_ld_aggregate.py` (orchestration with a fake reader) and `tests/genetics/test_ld_aggregate_cli.py` (subprocess)

The orchestration takes a `chunk_reader` callable so it can be unit-tested without Hail; the CLI passes `ld_hail.read_a_block_chunk`.

- [ ] **Step 1: Write the failing test (orchestration with a fake reader)**

```python
# tests/genetics/test_ld_aggregate.py
def _chrom_chunk_meta(row_start, row_stop, col_start, col_stop, chunk_dir):
    return {"row_start": row_start, "row_stop": row_stop,
            "column_start": col_start, "column_stop": col_stop,
            "directory": chunk_dir}


def test_accumulate_one_cohort_builds_precursor(tmp_path):
    cohort = tmp_path / "cohort"
    cohort.mkdir()
    # two variants on chr1 within 1Mb of each other
    rows = [("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0),
            ("1", 200, "C", "T", "1:200:C:T", 4, 0, 0.5)]
    import gzip
    with gzip.open(cohort / "variants.tsv.gz", "wt") as fh:
        fh.write("chr\tpos\tref\talt\tvariant_id\tn_nonmissing\tn_imputed\tgenotype_mean\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    np.save(cohort / "B.npy", np.array([[4.0, 90.0, 1.5], [2.0, 70.0, 1.0]]))
    np.save(cohort / "D.npy", np.array([[4.0, 180, 6], [180, 8200, 270], [6, 270, 10]]))

    manifest = {
        "study_name": "cohortA", "genome_build": "GRCh37",
        "covariate_schema": {"schema_id": "intercept_age_sex",
            "matrix_columns": ["intercept", "Age_numeric", "Sex_factor"],
            "sex_factor_recode": {"M": 1.0, "F": 2.0}},
        "variant_index": {"schema_version": "v0.3-with-genotype-stats"},
        "A_blocks": {"radius_bp": 1_000_000, "block_size": 8,
            "chromosomes": {"1": {"chunks": [
                _chrom_chunk_meta(0, 2, 0, 2, "A_blocks/chr1/chunk_000000")]}}},
    }
    import json
    (cohort / "manifest.json").write_text(json.dumps(manifest))

    # fake reader returns the dense A for the single chunk: X X^T for this cohort
    dense_a = np.array([[20.0, 8.0], [8.0, 10.0]])  # symmetric; diag 20,10 offdiag 8
    def fake_reader(chunk_dir):
        return dense_a

    precursor = tmp_path / "precursor"
    agg.accumulate(cohort, precursor, chunk_reader=fake_reader, pair_batch_rows=16)

    pm = agg.read_precursor_manifest(precursor)
    assert pm["n_cohorts"] == 1
    assert pm["cohorts"][0]["study_name"] == "cohortA"
    table = agg.read_variant_table(precursor).set_index("variant_id")
    assert table.loc["1:100:G:A", "a_diag"] == 20.0
    assert table.loc["1:200:C:T", "a_diag"] == 10.0
    import pyarrow.parquet as pq
    pairs = pq.read_table(precursor / "A_pairs" / "chr1.parquet").to_pandas()
    assert pairs["value"].tolist() == [8.0]
    np.testing.assert_array_equal(np.load(precursor / "D.npy"),
                                  np.load(cohort / "D.npy"))


def test_accumulate_rejects_duplicate_study(tmp_path):
    # (reuse the cohort built above via a helper; abbreviated here)
    ...  # build cohort as in the previous test
    # second accumulate of the same study_name must raise
    with pytest.raises(ValueError, match="already accumulated"):
        agg.accumulate(cohort, precursor, chunk_reader=fake_reader)
```

> Implementer note: factor the cohort-building into a `_build_synthetic_cohort(tmp_path)` helper in the test file and call it from both tests rather than duplicating.

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k accumulate -v`
Expected: FAIL with `AttributeError: ... 'accumulate'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py
import time
from datetime import datetime, timezone

LOCK_STALE_SECONDS = 24 * 3600


def _acquire_lock(precursor_dir: Path) -> Path:
    lock = precursor_paths(precursor_dir)["lock"]
    lock.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        os.write(fd, str(os.getpid()).encode())
        os.close(fd)
    except FileExistsError as exc:
        raise RuntimeError(
            f"Another accumulate is in progress (lock at {lock}). "
            "Remove it only if no process is running."
        ) from exc
    return lock


def _grouped_chrom_arrays(table: pd.DataFrame) -> dict[str, dict]:
    """Per-chromosome canonical-order positions + stable_ids for the new table."""
    out: dict[str, dict] = {}
    for chrom, sub in table.sort_values(["chr", "pos", "ref", "alt"]).groupby("chr", sort=False):
        out[chrom] = {
            "variant_id": sub["variant_id"].to_numpy(),
            "pos": sub["pos"].to_numpy(np.int64),
            "sid": sub["stable_id"].to_numpy(np.int64),
        }
    return out


def accumulate(
    cohort_dir: str | Path,
    precursor_dir: str | Path,
    chunk_reader,
    pair_batch_rows: int = 1_000_000,
    force: bool = False,
) -> None:
    """Add one cohort's A/B/D into the precursor (see spec accumulate phase)."""
    cohort_dir = Path(cohort_dir)
    precursor_dir = Path(precursor_dir)
    cohort_manifest = json.loads((cohort_dir / "manifest.json").read_text())
    study_name = cohort_manifest["study_name"]
    contract = extract_contract(cohort_manifest)

    pm = read_precursor_manifest(precursor_dir)
    if pm is None:
        pm = {
            "schema_version": PRECURSOR_SCHEMA_VERSION,
            "panel_specification_version": DEFAULT_PANEL_SPEC_VERSION,
            "contract": contract, "n_cohorts": 0, "cohorts": [],
            "next_stable_id": 0,
        }
        np.save(precursor_paths(precursor_dir)["d"], np.zeros((3, 3)))
    else:
        validate_contract(pm["contract"], contract)
        if any(c["study_name"] == study_name for c in pm["cohorts"]) and not force:
            raise ValueError(
                f"Study '{study_name}' already accumulated into this precursor; "
                "rebuild or pass force=True to override."
            )

    lock = _acquire_lock(precursor_dir)
    try:
        variants, b_mat, d_mat = read_cohort_assets(cohort_dir)
        table = read_variant_table(precursor_dir)
        table, id_map, next_id = merge_variant_table(
            table, variants, b_mat, pm["next_stable_id"])

        # cohort per-chrom canonical arrays for index -> (pos, sid) mapping
        cohort_chrom = {}
        for chrom, sub in variants.groupby("chr", sort=False):
            cohort_chrom[chrom] = {
                "pos": sub["pos"].to_numpy(np.int64),
                "sid": np.array([id_map[v] for v in sub["variant_id"]], dtype=np.int64),
            }

        diag_acc: dict[int, float] = {}
        for chrom, cmeta in cohort_manifest["A_blocks"]["chromosomes"].items():
            positions = cohort_chrom[chrom]["pos"]
            sids = cohort_chrom[chrom]["sid"]
            off_parts = []
            for chunk in cmeta["chunks"]:
                dense = chunk_reader(cohort_dir / chunk["directory"])
                diag_sids, diag_vals, offdiag = extract_chunk_entries(
                    dense, chunk["row_start"], chunk["column_start"],
                    positions, sids, contract["radius_bp"])
                for sid, val in zip(diag_sids.tolist(), diag_vals.tolist()):
                    diag_acc[sid] = diag_acc.get(sid, 0.0) + val
                if offdiag.size:
                    off_parts.append(offdiag)
            if off_parts:
                incoming = np.concatenate(off_parts)
                merge_pairs_into_file(
                    precursor_paths(precursor_dir)["pairs_dir"] / f"chr{chrom}.parquet",
                    incoming, batch_rows=pair_batch_rows)

        # apply diagonal accumulation
        if diag_acc:
            add = table["stable_id"].map(lambda s: diag_acc.get(int(s), 0.0))
            table["a_diag"] = table["a_diag"] + add
        write_variant_table(precursor_dir, table)

        d_path = precursor_paths(precursor_dir)["d"]
        np.save(d_path, np.load(d_path) + d_mat)

        pm["next_stable_id"] = next_id
        pm["n_cohorts"] += 1
        pm["cohorts"].append({
            "study_name": study_name,
            "accumulated_at_utc": datetime.now(timezone.utc).isoformat(),
            "n_variants": int(len(variants)),
        })
        write_precursor_manifest(precursor_dir, pm)  # commit point
    finally:
        lock.unlink(missing_ok=True)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k accumulate -v`
Expected: PASS.

- [ ] **Step 5: Write the CLI (rewrite `ld_aggregate_stats.py`, accumulate branch) and a subprocess test**

```python
# tests/genetics/test_ld_aggregate_cli.py
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from synthetic_plink import write_synthetic_plink

pytest.importorskip("hail", reason="hail not installed")
REPO_ROOT = Path(__file__).resolve().parents[2]
PREP = REPO_ROOT / "resources" / "genetics" / "ld_prepare_stats.py"
AGG = REPO_ROOT / "resources" / "genetics" / "ld_aggregate_stats.py"


def _run_prepare(cohort, calls, variants, samples, covs_rows, tmp):
    d = tmp / cohort
    d.mkdir()
    write_synthetic_plink(d / "geno", calls=calls, variants=variants, samples=samples)
    (d / "cov.txt").write_text("IID\tAge_numeric\tSex_factor\n" +
        "".join(f"{i}\t{a}\t{s}\n" for i, a, s in covs_rows))
    out = d / "out"
    subprocess.run([sys.executable, str(PREP), "--study-name", cohort,
        "--bfile", str(d / "geno"), "--covariates", str(d / "cov.txt"),
        "--output-dir", str(out), "--log-file", str(d / "log.txt")], check=True)
    return out


def test_cli_accumulate_creates_precursor(tmp_path):
    variants = [("22", 1_000_000, "A", "G"), ("22", 1_200_000, "C", "T")]
    out = _run_prepare("cohortA", [[0, 1, 2], [1, 0, 2]], variants,
        ["A1", "A2", "A3"], [("A1", 30, "M"), ("A2", 40, "F"), ("A3", 50, "M")], tmp_path)
    precursor = tmp_path / "precursor"
    r = subprocess.run([sys.executable, str(AGG), "--mode", "accumulate",
        "--cohort-dir", str(out), "--precursor-dir", str(precursor),
        "--log-file", str(tmp_path / "agg.log")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert (precursor / "precursor_manifest.json").is_file()
    assert json.loads((precursor / "precursor_manifest.json").read_text())["n_cohorts"] == 1
```

```python
# resources/genetics/ld_aggregate_stats.py  (full rewrite — accumulate branch)
#!/usr/bin/env python
import argparse
from pathlib import Path

import ld_aggregate as agg
from ld_hail import init_hail, read_a_block_chunk


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Section-15b central LD aggregation")
    p.add_argument("--mode", required=True, choices=["accumulate", "finalise"])
    p.add_argument("--precursor-dir", required=True)
    p.add_argument("--log-file", required=True)
    p.add_argument("--cohort-dir", help="15a cohort output dir (accumulate mode)")
    p.add_argument("--force", action="store_true")
    p.add_argument("--panel-dir", help="output panel dir (finalise mode)")
    p.add_argument("--maf-threshold", type=float, default=agg.DEFAULT_MAF_THRESHOLD)
    p.add_argument("--min-adj-diag", type=float, default=agg.DEFAULT_MIN_ADJ_DIAG)
    p.add_argument("--min-cohorts", type=int, default=None)
    p.add_argument("--a-block-size", type=int, default=None)
    p.add_argument("--a-max-dense-gb", type=float, default=1.0)
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "accumulate":
        if not args.cohort_dir:
            raise SystemExit("--cohort-dir is required for accumulate mode")
        hail_log = Path(args.log_file).parent / "hail.log"
        init_hail(hail_log)
        agg.accumulate(args.cohort_dir, args.precursor_dir,
                       chunk_reader=read_a_block_chunk, force=args.force)
        print(f"[section15b] accumulated cohort into {args.precursor_dir}", flush=True)
    else:
        # finalise branch added in Task 11
        raise SystemExit("finalise mode not yet implemented")


if __name__ == "__main__":
    main()
```

- [ ] **Step 6: Run the CLI test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate_cli.py -k accumulate -v`
Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add resources/genetics/ld_aggregate.py resources/genetics/ld_aggregate_stats.py tests/genetics/test_ld_aggregate.py tests/genetics/test_ld_aggregate_cli.py
git commit -m "feat(15b): accumulate orchestration + CLI accumulate mode"
```

---

## Task 9: Intersection resolution + pooled filters

**Files:**
- Modify: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

`resolve_intersection` keeps variants with `membership_count >= min_cohorts` (default `K`), sorts by `(chr,pos,ref,alt)`, assigns contiguous `pooled_index`. `apply_filters` computes pooled MAF and the adjusted diagonal, drops failing variants, returns the kept table plus a dropped-reasons frame.

- [ ] **Step 1: Write the failing test**

```python
def test_resolve_intersection_keeps_full_membership():
    table = pd.DataFrame({
        "stable_id": [0, 1, 2], "variant_id": ["1:200:C:T", "1:100:G:A", "2:50:A:C"],
        "chr": ["1", "1", "2"], "pos": [200, 100, 50], "ref": ["C", "G", "A"],
        "alt": ["T", "A", "C"], "membership_count": [2, 1, 2],
        "b_intercept": [0.0]*3, "b_age": [0.0]*3, "b_sex": [0.0]*3,
        "a_diag": [0.0]*3, "n_nonmissing": [0]*3, "n_imputed": [0]*3,
    })
    kept = agg.resolve_intersection(table, n_cohorts=2, min_cohorts=2)
    # variant 1 (membership 1) dropped; survivors sorted by chr,pos
    assert kept["variant_id"].tolist() == ["1:200:C:T", "2:50:A:C"]
    assert kept["pooled_index"].tolist() == [0, 1]


def test_apply_filters_drops_rare_and_unstable():
    # D[0,0] = N = 100. b_intercept 30 -> f=0.15 keep; b_intercept 1 -> f=0.005 drop
    d = np.diag([100.0, 1.0, 1.0])
    kept = pd.DataFrame({
        "stable_id": [0, 1], "variant_id": ["1:100:G:A", "1:200:C:T"],
        "chr": ["1", "1"], "pos": [100, 200], "ref": ["G", "C"], "alt": ["A", "T"],
        "b_intercept": [30.0, 1.0], "b_age": [0.0, 0.0], "b_sex": [0.0, 0.0],
        "a_diag": [50.0, 50.0], "pooled_index": [0, 1],
    })
    surv, dropped = agg.apply_filters(kept, d, maf_threshold=0.01, min_adj_diag=0.0)
    assert surv["variant_id"].tolist() == ["1:100:G:A"]
    assert surv["pooled_index"].tolist() == [0]                 # re-indexed
    assert "1:200:C:T" in dropped["variant_id"].tolist()
    assert dropped.set_index("variant_id").loc["1:200:C:T", "reason"] == "maf_below_threshold"
    assert "a_adj_diag" in surv.columns
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k "intersection or apply_filters" -v`
Expected: FAIL with `AttributeError: ... 'resolve_intersection'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py

def resolve_intersection(table: pd.DataFrame, n_cohorts: int,
                         min_cohorts: int | None = None) -> pd.DataFrame:
    """Keep variants present in enough cohorts; assign contiguous pooled_index."""
    threshold = n_cohorts if min_cohorts is None else min_cohorts
    kept = table[table["membership_count"] >= threshold].copy()
    kept = kept.sort_values(["chr", "pos", "ref", "alt"],
                            key=lambda s: s.map(int) if s.name == "chr" else s)
    kept = kept.reset_index(drop=True)
    kept["pooled_index"] = np.arange(len(kept), dtype=np.int64)
    return kept


def apply_filters(kept: pd.DataFrame, d_matrix: np.ndarray,
                  maf_threshold: float, min_adj_diag: float
                  ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Drop rare variants and unstable adjusted diagonals; re-index survivors.

    Returns (survivors with pooled_index reset + a_adj_diag column,
    dropped frame with variant_id + reason).
    """
    n = d_matrix[0, 0]
    d_inv = np.linalg.inv(d_matrix)
    b = kept[["b_intercept", "b_age", "b_sex"]].to_numpy()
    freq = b[:, 0] / (2.0 * n)
    maf = np.minimum(freq, 1.0 - freq)
    # adjusted diagonal a_diag - b D^-1 b^T per row
    adj_diag = kept["a_diag"].to_numpy() - np.einsum("ij,jk,ik->i", b, d_inv, b)

    reason = np.where(maf < maf_threshold, "maf_below_threshold",
              np.where(adj_diag <= min_adj_diag, "unstable_adj_diagonal", ""))
    drop_mask = reason != ""
    dropped = pd.DataFrame({
        "variant_id": kept.loc[drop_mask, "variant_id"].to_numpy(),
        "reason": reason[drop_mask],
    })
    surv = kept.loc[~drop_mask].copy()
    surv["a_adj_diag"] = adj_diag[~drop_mask]
    surv = surv.sort_values("pooled_index").reset_index(drop=True)
    surv["pooled_index"] = np.arange(len(surv), dtype=np.int64)
    return surv, dropped
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k "intersection or apply_filters" -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): intersection resolution and pooled MAF/diagonal filters"
```

---

## Task 10: Per-chunk A_adj and R block math

**Files:**
- Modify: `resources/genetics/ld_aggregate.py`
- Test: `tests/genetics/test_ld_aggregate.py`

`compute_r_block` takes a dense pooled `A` window block (off-diagonal scattered + diagonal placed), the pooled `B`/`W` slices, and the `a_adj_diag` slices, and returns the `R` block plus the per-row sparsify intervals. Verified against a direct full-matrix residualisation on a tiny example.

- [ ] **Step 1: Write the failing test**

```python
def test_compute_r_block_matches_direct_residualisation():
    rng = np.random.default_rng(0)
    n, p = 40, 5
    X = rng.integers(0, 3, size=(n, p)).astype(float)
    C = np.column_stack([np.ones(n), rng.normal(size=n), rng.integers(1, 3, n)])
    A = X.T @ X
    B = X.T @ C
    D = C.T @ C
    d_inv = np.linalg.inv(D)
    A_adj = A - B @ d_inv @ B.T
    diag = np.diag(A_adj)
    R_ref = A_adj / np.sqrt(np.outer(diag, diag))

    # full window: rows 0..p, cols 0..p, all pairs in-window
    W = B @ d_inv
    a_block = A.copy()                       # dense pooled A window (incl diagonal)
    starts = np.arange(p, dtype=np.int64)    # upper-tri row intervals
    stops = np.full(p, p, dtype=np.int64)
    r_block, out_starts, out_stops = agg.compute_r_block(
        a_block, B, W, diag, diag, row_start=0, col_start=0,
        starts=starts, stops=stops)
    np.testing.assert_allclose(np.triu(r_block), np.triu(R_ref), rtol=1e-10)
    assert out_starts.tolist() == starts.tolist()
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k compute_r_block -v`
Expected: FAIL with `AttributeError: ... 'compute_r_block'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py

def compute_r_block(
    a_block: np.ndarray,
    b_rows: np.ndarray,
    w_cols: np.ndarray,
    adj_diag_rows: np.ndarray,
    adj_diag_cols: np.ndarray,
    row_start: int,
    col_start: int,
    starts: np.ndarray,
    stops: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert a dense pooled A window block into an R block.

    a_block[r, c] is pooled A for pooled row row_start+r, pooled col col_start+c.
    b_rows: B for the block's rows (n_rows x 3). w_cols: W=(B D^-1) for the block's
    columns (n_cols x 3). adj_diag_* are the A_adj diagonals for rows/cols.
    starts/stops are chunk-local column intervals (upper triangle within window).
    """
    corr = b_rows @ w_cols.T                       # rank-3 covariate correction
    a_adj = a_block - corr
    denom = np.sqrt(np.outer(adj_diag_rows, adj_diag_cols))
    r_block = a_adj / denom
    return r_block, np.asarray(starts, np.int64), np.asarray(stops, np.int64)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k compute_r_block -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add resources/genetics/ld_aggregate.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): per-chunk A_adj and R block math"
```

---

## Task 11: `finalise` orchestration — streaming, R_blocks, panel artefacts

**Files:**
- Modify: `resources/genetics/ld_aggregate.py` (finalise orchestration with injected R-writer)
- Modify: `resources/genetics/ld_aggregate_stats.py` (finalise branch)
- Test: `tests/genetics/test_ld_aggregate.py` (orchestration with a capture writer)

Finalise reuses bounded row-chunk planning. Add small pure helpers `_dense_product_gb` and `_bounded_row_stop` to `ld_aggregate.py` (same logic as `ld_hail`; duplicated deliberately to keep this module Hail-free). The R-writer is injected so the orchestration is unit-testable; the CLI passes `ld_hail.write_r_block_chunk`.

- [ ] **Step 1: Write the failing test (orchestration with a capture writer)**

```python
def test_finalise_writes_panel_and_calls_writer(tmp_path):
    # Build a 2-cohort precursor synthetically by accumulating two fake cohorts
    # via the Task 8 helper, then finalise with a capturing R-writer.
    precursor = _build_two_cohort_precursor(tmp_path)   # test helper (see note)
    panel = tmp_path / "panel"
    written = []
    def capture_writer(dense, starts, stops, out_dir, block_size):
        written.append((np.asarray(dense).copy(), Path(out_dir)))
        Path(out_dir).mkdir(parents=True, exist_ok=True)
        (Path(out_dir) / "MARKER").write_text("ok")

    agg.finalise(precursor, panel, r_writer=capture_writer,
                 maf_threshold=0.0, min_adj_diag=-1e9, block_size=8,
                 max_dense_gb=1.0)

    assert (panel / "pooled_manifest.json").is_file()
    assert (panel / "variants.tsv.gz").is_file()
    assert (panel / "dropped_variants.tsv").is_file()
    assert written, "R writer was never called"
    # diagonal of every emitted R block is 1 where present
    for dense, _ in written:
        d = np.diag(dense)
        np.testing.assert_allclose(d[d != 0], 1.0, rtol=1e-9)
```

> Implementer note: add `_build_two_cohort_precursor(tmp_path)` to the test file: build two synthetic cohorts (overlapping variant sets on one chromosome), `agg.accumulate` both with a fake `chunk_reader` returning each cohort's `X X^T`, and return the precursor path. This same helper anchors Task 12's e2e equivalence assertions.

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k finalise -v`
Expected: FAIL with `AttributeError: ... 'finalise'`.

- [ ] **Step 3: Write minimal implementation**

```python
# Append to ld_aggregate.py
import gzip as _gzip


def _dense_product_gb(n_rows: int, n_cols: int) -> float:
    return (n_rows * n_cols * 8) / (1024 ** 3)


def _bounded_row_stop(row_start, requested_row_stop, stops, max_dense_gb):
    best, low, high = None, row_start + 1, requested_row_stop
    while low <= high:
        mid = (low + high) // 2
        col_stop = int(stops[row_start:mid].max())
        if _dense_product_gb(mid - row_start, col_stop - row_start) <= max_dense_gb:
            best, low = mid, mid + 1
        else:
            high = mid - 1
    if best is None:
        raise MemoryError(
            f"finalise chunk exceeds max_dense_gb even for one row at {row_start}")
    return best


def finalise(precursor_dir, panel_dir, r_writer, maf_threshold, min_adj_diag,
             block_size, max_dense_gb=1.0, min_cohorts=None):
    """Resolve intersection, adjust, convert to R, write a versioned panel."""
    precursor_dir, panel_dir = Path(precursor_dir), Path(panel_dir)
    pm = read_precursor_manifest(precursor_dir)
    if pm is None or pm["n_cohorts"] == 0:
        raise ValueError("Cannot finalise an empty precursor")
    d_matrix = np.load(precursor_paths(precursor_dir)["d"])
    d_rank = int(np.linalg.matrix_rank(d_matrix))
    if d_rank < d_matrix.shape[0]:
        raise ValueError(f"Pooled D is rank-deficient (rank {d_rank}); cannot solve")

    table = read_variant_table(precursor_dir)
    kept = resolve_intersection(table, pm["n_cohorts"], min_cohorts)
    if len(kept) == 0:
        raise ValueError("Intersection is empty; no variant is present in all cohorts")
    surv, dropped = apply_filters(kept, d_matrix, maf_threshold, min_adj_diag)
    if len(surv) == 0:
        raise ValueError("All variants dropped by pooled filters")

    d_inv = np.linalg.inv(d_matrix)
    b_all = surv[["b_intercept", "b_age", "b_sex"]].to_numpy()
    w_all = b_all @ d_inv
    adj_diag_all = surv["a_adj_diag"].to_numpy()
    sid_to_pooled = dict(zip(surv["stable_id"].to_numpy(), surv["pooled_index"].to_numpy()))
    pos_all = surv["pos"].to_numpy(np.int64)

    panel_dir.mkdir(parents=True, exist_ok=True)
    r_root = panel_dir / "R_blocks"

    for chrom, sub in surv.groupby("chr", sort=False):
        idx = sub["pooled_index"].to_numpy(np.int64)
        local_pos = pos_all[idx]
        # local stops within this chromosome's pooled rows
        n_chr = len(idx)
        local_stops = np.searchsorted(local_pos, local_pos + pm["contract"]["radius_bp"],
                                      side="right")
        pair_path = precursor_paths(precursor_dir)["pairs_dir"] / f"chr{chrom}.parquet"
        pair_iter = _PairCursor(pair_path, sid_to_pooled, idx[0])

        row_start = 0
        chunk_no = 0
        while row_start < n_chr:
            requested = min(row_start + 50_000, n_chr)
            row_stop = _bounded_row_stop(row_start, requested, local_stops, max_dense_gb)
            col_stop = int(local_stops[row_start:row_stop].max())
            n_rows, n_cols = row_stop - row_start, col_stop - row_start
            block = np.zeros((n_rows, n_cols), dtype=np.float64)
            # diagonal
            for r in range(n_rows):
                block[r, (row_start + r) - row_start] = \
                    surv["a_diag"].to_numpy()[idx[row_start + r]]
            # off-diagonal from the sorted pair cursor, by chromosome-local pooled index
            pair_iter.fill_block(block, row_start, row_stop, col_stop, idx)

            g_rows = idx[row_start:row_stop]
            g_cols = idx[row_start:col_stop]
            starts = np.arange(n_rows, dtype=np.int64)               # local diag start
            stops = (local_stops[row_start:row_stop] - row_start).astype(np.int64)
            r_block, st, sp = compute_r_block(
                block, b_all[g_rows], w_all[g_cols],
                adj_diag_all[g_rows], adj_diag_all[g_cols],
                row_start, row_start, starts, stops)
            out_dir = r_root / f"chr{chrom}" / f"chunk_{chunk_no:06d}"
            r_writer(r_block, st, sp, out_dir, block_size)
            row_start, chunk_no = row_stop, chunk_no + 1

    _write_panel_artefacts(panel_dir, pm, surv, dropped, d_matrix, d_rank,
                           maf_threshold, min_adj_diag)


def _write_panel_artefacts(panel_dir, pm, surv, dropped, d_matrix, d_rank,
                           maf_threshold, min_adj_diag):
    cols = ["chr", "pos", "ref", "alt", "variant_id", "pooled_index",
            "b_intercept", "a_diag", "a_adj_diag"]
    surv["maf"] = np.minimum(surv["b_intercept"] / (2 * d_matrix[0, 0]),
                             1 - surv["b_intercept"] / (2 * d_matrix[0, 0]))
    with _gzip.open(panel_dir / "variants.tsv.gz", "wt") as fh:
        surv[cols + ["maf"]].to_csv(fh, sep="\t", index=False)
    dropped.to_csv(panel_dir / "dropped_variants.tsv", sep="\t", index=False)
    pd.DataFrame(pm["cohorts"]).to_csv(panel_dir / "cohort_inclusion.tsv",
                                       sep="\t", index=False)
    manifest = {
        "module": "15b", "panel_specification_version": pm["panel_specification_version"],
        "covariate_schema": pm["contract"]["schema_id"], "n_cohorts": pm["n_cohorts"],
        "cohorts": [c["study_name"] for c in pm["cohorts"]],
        "n_variants_panel": int(len(surv)), "n_variants_dropped": int(len(dropped)),
        "d_rank": d_rank, "d_condition_number": float(np.linalg.cond(d_matrix)),
        "maf_threshold": maf_threshold, "min_adj_diag": min_adj_diag,
        "regularisation": "none", "r_dtype": "float32",
    }
    (panel_dir / "pooled_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (panel_dir / "qc_report.txt").write_text(
        f"Section 15b panel built: {len(surv)} variants, {len(dropped)} dropped, "
        f"{pm['n_cohorts']} cohorts, D rank {d_rank}.\n")
```

```python
# Append to ld_aggregate.py — streaming cursor over a sorted pair file

class _PairCursor:
    """Sequential reader over a position-sorted A_pairs/chr*.parquet file.

    fill_block scatters off-diagonal pooled entries whose row pooled-index is in
    the current chunk into the dense block, mapping sid -> pooled index. Pairs
    with an endpoint outside the intersection (sid not in the map) are skipped.
    """
    def __init__(self, path, sid_to_pooled, chrom_first_pooled):
        self._sid_to_pooled = sid_to_pooled
        self._rows = (pq.read_table(path).to_pandas()
                      if Path(path).is_file() else None)
        self._cursor = 0

    def fill_block(self, block, row_start, row_stop, col_stop, idx):
        if self._rows is None:
            return
        n = len(self._rows)
        pos_i = self._rows["pos_i"].to_numpy()
        sid_i = self._rows["sid_i"].to_numpy()
        sid_j = self._rows["sid_j"].to_numpy()
        val = self._rows["value"].to_numpy()
        # pooled row indices for this chromosome chunk are idx[row_start:row_stop]
        pooled_rows = set(range(row_start, row_stop))
        while self._cursor < n:
            pi = self._sid_to_pooled.get(int(sid_i[self._cursor]))
            pj = self._sid_to_pooled.get(int(sid_j[self._cursor]))
            if pi is None or pj is None:
                self._cursor += 1
                continue
            if pi >= row_stop:
                break                                  # past this chunk's rows
            if pi < row_start:
                self._cursor += 1
                continue
            block[pi - row_start, pj - row_start] += val[self._cursor]
            self._cursor += 1
```

> Implementer note: the `_PairCursor` relies on chromosome-global pooled indices throughout, and on the pair file being sorted by `pos_i` so the pooled row index is non-decreasing. Verify the cursor advances monotonically (does not revisit rows across chunks) with the Task 11 test before wiring the CLI.

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate.py -k finalise -v`
Expected: PASS.

- [ ] **Step 5: Wire the CLI finalise branch**

```python
# In ld_aggregate_stats.py main(), replace the finalise branch:
    else:
        if not args.panel_dir:
            raise SystemExit("--panel-dir is required for finalise mode")
        hail_log = Path(args.log_file).parent / "hail.log"
        init_hail(hail_log)
        from ld_hail import write_r_block_chunk, DEFAULT_A_BLOCK_SIZE
        agg.finalise(
            args.precursor_dir, args.panel_dir, r_writer=write_r_block_chunk,
            maf_threshold=args.maf_threshold, min_adj_diag=args.min_adj_diag,
            block_size=args.a_block_size or DEFAULT_A_BLOCK_SIZE,
            max_dense_gb=args.a_max_dense_gb, min_cohorts=args.min_cohorts)
        print(f"[section15b] finalised panel into {args.panel_dir}", flush=True)
```

- [ ] **Step 6: Commit**

```bash
git add resources/genetics/ld_aggregate.py resources/genetics/ld_aggregate_stats.py tests/genetics/test_ld_aggregate.py
git commit -m "feat(15b): finalise orchestration, R_blocks streaming, panel artefacts"
```

---

## Task 12: End-to-end synthetic two-cohort equivalence

**Files:**
- Modify: `tests/genetics/test_ld_aggregate_synthetic.py`
- Test: same file (drop the skip + add equivalence/order/guard tests)

This is the statistical acceptance gate. Run two synthetic cohorts through real 15a (`ld_prepare_stats.py`), `accumulate` both via the real CLI (real Hail `read_a_block_chunk`), `finalise`, read the `R_blocks` back, and assert equality with a direct stacked-residualisation NumPy reference on the intersection.

- [ ] **Step 1: Replace the skipped test with the real assertion**

```python
# tests/genetics/test_ld_aggregate_synthetic.py — replace the skip body
import subprocess, sys
from pathlib import Path
import numpy as np
from hail.linalg import BlockMatrix

AGG = REPO_ROOT / "resources" / "genetics" / "ld_aggregate_stats.py"


def _accumulate(cohort_out, precursor, tmp_path):
    r = subprocess.run([sys.executable, str(AGG), "--mode", "accumulate",
        "--cohort-dir", str(cohort_out), "--precursor-dir", str(precursor),
        "--log-file", str(tmp_path / "agg.log")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def _finalise(precursor, panel, tmp_path):
    r = subprocess.run([sys.executable, str(AGG), "--mode", "finalise",
        "--precursor-dir", str(precursor), "--panel-dir", str(panel),
        "--maf-threshold", "0.0", "--min-adj-diag", "-1e9",
        "--a-block-size", "8", "--log-file", str(tmp_path / "fin.log")],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def _read_r_full(panel, chrom, n):
    """Reassemble the full pooled R matrix for a chromosome from its chunks."""
    import json
    out = np.zeros((n, n))
    chunk_root = panel / "R_blocks" / f"chr{chrom}"
    starts = sorted(int(p.name.split("_")[1]) for p in chunk_root.glob("chunk_*"))
    row = 0
    for c in sorted(chunk_root.glob("chunk_*")):
        block = BlockMatrix.read(str(c)).to_numpy()
        out[row:row + block.shape[0], row:row + block.shape[1]] += block
        row += block.shape[0]
    return out


def test_two_cohort_a_adj_matches_stacked_residualisation(tmp_path):
    out_a = _run_15a(cohort_name="cohort_a", calls=COHORT_A_CALLS,
        variants=SHARED_VARIANTS, samples=COHORT_A_SAMPLES,
        covariates_rows=COHORT_A_COVARIATES, tmp_path=tmp_path)
    out_b = _run_15a(cohort_name="cohort_b", calls=COHORT_B_CALLS,
        variants=SHARED_VARIANTS, samples=COHORT_B_SAMPLES,
        covariates_rows=COHORT_B_COVARIATES, tmp_path=tmp_path)

    precursor, panel = tmp_path / "precursor", tmp_path / "panel"
    _accumulate(out_a, precursor, tmp_path)
    _accumulate(out_b, precursor, tmp_path)
    _finalise(precursor, panel, tmp_path)

    # direct reference on the shared variant set (all SHARED_VARIANTS overlap)
    x_full = np.vstack([np.asarray(COHORT_A_CALLS, float).T,
                        np.asarray(COHORT_B_CALLS, float).T])
    c_full = np.vstack([_covariate_matrix(COHORT_A_COVARIATES),
                        _covariate_matrix(COHORT_B_COVARIATES)])
    A, B, D = x_full.T @ x_full, x_full.T @ c_full, c_full.T @ c_full
    A_adj = A - B @ np.linalg.inv(D) @ B.T
    diag = np.diag(A_adj)
    R_ref = A_adj / np.sqrt(np.outer(diag, diag))

    R_got = _read_r_full(panel, "22", len(SHARED_VARIANTS))
    np.testing.assert_allclose(np.triu(R_got), np.triu(R_ref), rtol=1e-6, atol=1e-6)
```

- [ ] **Step 2: Run it to verify it passes**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate_synthetic.py -v`
Expected: PASS (both the existing D/B test and the new A_adj/R test).

- [ ] **Step 3: Add order-independence, intersection, and guard tests**

```python
def test_accumulate_order_independent(tmp_path):
    out_a = _run_15a(cohort_name="cohort_a", calls=COHORT_A_CALLS,
        variants=SHARED_VARIANTS, samples=COHORT_A_SAMPLES,
        covariates_rows=COHORT_A_COVARIATES, tmp_path=tmp_path)
    out_b = _run_15a(cohort_name="cohort_b", calls=COHORT_B_CALLS,
        variants=SHARED_VARIANTS, samples=COHORT_B_SAMPLES,
        covariates_rows=COHORT_B_COVARIATES, tmp_path=tmp_path)
    p_ab, p_ba = tmp_path / "p_ab", tmp_path / "p_ba"
    _accumulate(out_a, p_ab, tmp_path); _accumulate(out_b, p_ab, tmp_path)
    _accumulate(out_b, p_ba, tmp_path); _accumulate(out_a, p_ba, tmp_path)
    panel_ab, panel_ba = tmp_path / "ab", tmp_path / "ba"
    _finalise(p_ab, panel_ab, tmp_path); _finalise(p_ba, panel_ba, tmp_path)
    np.testing.assert_allclose(_read_r_full(panel_ab, "22", 3),
                               _read_r_full(panel_ba, "22", 3), rtol=1e-9)


def test_intersection_excludes_partial_variant(tmp_path):
    # cohort B is missing the middle variant -> it must not appear in the panel
    out_a = _run_15a(cohort_name="cohort_a", calls=COHORT_A_CALLS,
        variants=SHARED_VARIANTS, samples=COHORT_A_SAMPLES,
        covariates_rows=COHORT_A_COVARIATES, tmp_path=tmp_path)
    partial = [SHARED_VARIANTS[0], SHARED_VARIANTS[2]]
    calls_b = [COHORT_B_CALLS[0], COHORT_B_CALLS[2]]
    out_b = _run_15a(cohort_name="cohort_b", calls=calls_b, variants=partial,
        samples=COHORT_B_SAMPLES, covariates_rows=COHORT_B_COVARIATES, tmp_path=tmp_path)
    precursor, panel = tmp_path / "pp", tmp_path / "panel_partial"
    _accumulate(out_a, precursor, tmp_path); _accumulate(out_b, precursor, tmp_path)
    _finalise(precursor, panel, tmp_path)
    import gzip
    with gzip.open(panel / "variants.tsv.gz", "rt") as fh:
        body = fh.read()
    assert "2000000" not in body and "3000000" in body and "1000000" in body


def test_accumulate_duplicate_study_refused(tmp_path):
    out_a = _run_15a(cohort_name="cohort_a", calls=COHORT_A_CALLS,
        variants=SHARED_VARIANTS, samples=COHORT_A_SAMPLES,
        covariates_rows=COHORT_A_COVARIATES, tmp_path=tmp_path)
    precursor = tmp_path / "dup"
    _accumulate(out_a, precursor, tmp_path)
    r = subprocess.run([sys.executable, str(AGG), "--mode", "accumulate",
        "--cohort-dir", str(out_a), "--precursor-dir", str(precursor),
        "--log-file", str(tmp_path / "agg2.log")], capture_output=True, text=True)
    assert r.returncode != 0
    assert "already accumulated" in (r.stderr + r.stdout)
```

- [ ] **Step 4: Run the full synthetic suite**

Run: `conda run -n hail_env python -m pytest tests/genetics/test_ld_aggregate_synthetic.py -v`
Expected: PASS (all tests; no skips).

- [ ] **Step 5: Commit**

```bash
git add tests/genetics/test_ld_aggregate_synthetic.py
git commit -m "test(15b): end-to-end two-cohort A_adj/R equivalence, order, intersection, guards"
```

---

## Task 13: Pipeline integration — stage script, parameters, checks, wiki

**Files:**
- Rewrite: `15b-ld_aggregate_stats.sh`
- Modify: `resources/parameters`
- Modify: `resources/logs/check_logs.sh`, `resources/logs/check_results.sh`
- Modify: `godmc_phase2.wiki/Run-federated-LD-reference-panel.md`

- [ ] **Step 1: Add parameters**

In `resources/parameters`, after the existing `ld_*` block (around line 256), add:

```bash
ld_precursor_dir="${ld_dir}/precursor"
ld_panel_dir="${ld_dir}/panels"
ld_15b_mode="${ld_15b_mode:-accumulate}"
ld_maf_threshold="${ld_maf_threshold:-0.01}"
ld_min_adj_diag="${ld_min_adj_diag:-0.0}"
ld_min_cohorts="${ld_min_cohorts:-}"
```

- [ ] **Step 2: Rewrite the stage script**

```bash
# 15b-ld_aggregate_stats.sh
#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_15b_logfile})
print_version

mkdir -p ${ld_dir}
mkdir -p ${ld_precursor_dir}
mkdir -p ${ld_panel_dir}

if [ "${ld_15b_mode}" = "accumulate" ]; then
	if [ ! -f "${ld_prepare_dir}/manifest.json" ]; then
		echo "Problem: cohort 15a manifest is required at ${ld_prepare_dir}/manifest.json"
		exit 1
	fi
	python resources/genetics/ld_aggregate_stats.py \
		--mode accumulate \
		--cohort-dir ${ld_prepare_dir} \
		--precursor-dir ${ld_precursor_dir} \
		--log-file ${section_15b_logfile}
	echo "Successfully accumulated cohort into the LD precursor"
elif [ "${ld_15b_mode}" = "finalise" ]; then
	min_cohorts_arg=""
	if [ -n "${ld_min_cohorts}" ]; then
		min_cohorts_arg="--min-cohorts ${ld_min_cohorts}"
	fi
	python resources/genetics/ld_aggregate_stats.py \
		--mode finalise \
		--precursor-dir ${ld_precursor_dir} \
		--panel-dir ${ld_panel_dir} \
		--maf-threshold ${ld_maf_threshold} \
		--min-adj-diag ${ld_min_adj_diag} \
		${min_cohorts_arg} \
		--log-file ${section_15b_logfile}
	echo "Successfully finalised the pooled LD panel"
else
	echo "Problem: ld_15b_mode must be 'accumulate' or 'finalise'; got '${ld_15b_mode}'"
	exit 1
fi
```

- [ ] **Step 3: Update log/result checks**

In `resources/logs/check_logs.sh`, update the 15b success-marker check to accept either completion line:

```bash
# 15b: accept accumulate or finalise success
grep -E "Successfully (accumulated cohort into the LD precursor|finalised the pooled LD panel)" \
	${section_15b_logfile} > /dev/null
```

In `resources/logs/check_results.sh`, extend `check_results_15` to assert the precursor manifest after accumulate and the panel manifest after finalise:

```bash
# after the existing cohort-side asserts in check_results_15
if [ -d "${ld_precursor_dir}" ]; then
	check_file_exists "${ld_precursor_dir}/precursor_manifest.json"
fi
if [ -d "${ld_panel_dir}" ] && [ -n "$(ls -A ${ld_panel_dir} 2>/dev/null)" ]; then
	check_file_exists "${ld_panel_dir}/pooled_manifest.json"
fi
```

> Implementer note: read the existing `check_results_15` body first and match its helper names (`check_file_exists` or the project's equivalent) exactly.

- [ ] **Step 4: Update the wiki**

In `godmc_phase2.wiki/Run-federated-LD-reference-panel.md`, add a "Central aggregation (15b)" section documenting the two-phase workflow:

```markdown
## Central aggregation (section 15b)

15b runs centrally in two modes, selected by `ld_15b_mode`:

- **accumulate** (default): add one cohort's 15a output (`${ld_prepare_dir}`) into the
  running precursor at `${ld_precursor_dir}`. Run once per cohort as cohorts arrive; a
  cohort's raw `A_blocks` can be deleted after it is accumulated. Re-accumulating the same
  `study_name` is refused.
- **finalise**: build an immutable pooled LD panel under `${ld_panel_dir}` from the current
  precursor — intersection of variants present in all cohorts, pooled MAF filter
  (`ld_maf_threshold`, default 0.01), covariate adjustment, and float32 `R_blocks`.

Example:

    ld_15b_mode=accumulate ./15b-ld_aggregate_stats.sh <config>   # per cohort
    ld_15b_mode=finalise   ./15b-ld_aggregate_stats.sh <config>   # once cohorts are in
```

- [ ] **Step 5: Verify the scaffold is fully replaced**

Run: `conda run -n hail_env python -m pytest tests/genetics/ -v`
Expected: PASS (whole section-15 suite, no skips).

Also confirm no caller references the removed scaffold behaviour:

Run: `grep -rn "scaffold" resources/genetics/ld_aggregate_stats.py 15b-ld_aggregate_stats.sh`
Expected: no matches.

- [ ] **Step 6: Commit**

```bash
git add 15b-ld_aggregate_stats.sh resources/parameters resources/logs/check_logs.sh resources/logs/check_results.sh godmc_phase2.wiki/Run-federated-LD-reference-panel.md
git commit -m "feat(15b): wire two-phase aggregation into the pipeline + docs"
```

---

## Post-implementation: update the progress tracker

After Task 13, update `LD_IMPLEMENTATION_PROGRESS.md`: move the 15b milestone rows (manifest validation, pooled variant index, block alignment, aggregation, covariate adjustment, final R, numerical policy, QC reporting) to Done/In-progress as appropriate, record the decisions taken in this plan (two-phase incremental precursor, position-sorted pair store, intersection semantics, no checksums, MAF default 0.01, float32 R), and note the two-cohort 15b pilot as the next scale step. Commit separately.

---

## Self-review notes (for the executing agent)

- **Spec coverage:** accumulate phase → Tasks 4–8; finalise phase → Tasks 9–11; intersection → Task 9; A_adj/R math → Task 10; error handling/guards → Tasks 2, 8, 11, 12; validation tests → Task 12; pipeline integration → Task 13. The pair-store streaming refinement (sorted by genomic position) is realised in Tasks 5–6 and consumed in Task 11's `_PairCursor`.
- **Out of scope (per spec):** globally projected PCs, cohort fixed effects, global ridge, batching the *incoming* pair side for genome-wide scale (release-1 note in Task 6).
- **Test helpers to factor (not duplicate):** `_build_synthetic_cohort` (Task 8), `_build_two_cohort_precursor` (Task 11) — define once in the test file and reuse; Task 8's duplicate-study test reuses the former, Task 11's finalise test reuses the latter.
</content>
