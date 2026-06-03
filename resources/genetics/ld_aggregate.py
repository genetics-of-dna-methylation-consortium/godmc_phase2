# resources/genetics/ld_aggregate.py
"""Pure-Python central LD aggregation for section 15b.

This module deliberately does NOT import hail. Every code path that touches
Hail BlockMatrix IO lives in ld_hail.py; this module operates on NumPy arrays,
pandas DataFrames, and parquet files so its unit tests run without a JVM.
"""

from __future__ import annotations

import json
import os
from datetime import datetime, timezone
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


def _atomic_np_save(path: str | Path, arr: np.ndarray) -> None:
    """Atomically write a .npy file (temp-then-rename), like the other writers."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".tmp.{os.getpid()}.npy")
    np.save(tmp, arr)
    os.replace(tmp, path)


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
        precursor_dir.mkdir(parents=True, exist_ok=True)
        _atomic_np_save(precursor_paths(precursor_dir)["d"], np.zeros((3, 3)))
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
        _atomic_np_save(d_path, np.load(d_path) + d_mat)

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
