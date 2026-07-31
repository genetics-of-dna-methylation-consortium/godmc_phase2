# resources/genetics/ld_aggregate.py
"""Pure-Python central LD aggregation for section 15.

This module deliberately does NOT import hail. Every code path that touches
Hail BlockMatrix IO lives in ld_hail.py; this module operates on NumPy arrays,
pandas DataFrames, and parquet files so its unit tests run without a JVM.
"""

from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import ld_checksums

PRECURSOR_SCHEMA_VERSION = "15-aggregate-precursor-v1"
COHORT_VARIANT_STATS_SCHEMA_VERSION = "15-cohort-variant-stats-v1"
DEFAULT_PANEL_SPEC_VERSION = "0.1.0"
DEFAULT_MAF_THRESHOLD = 0.01
DEFAULT_MIN_ADJ_DIAG = 0.0

VARIANT_COLUMNS = [
    "stable_id", "variant_id", "chr", "pos", "ref", "alt",
    "membership_count", "b_intercept",
    "a_diag", "n_nonmissing", "n_imputed",
]

_VARIANT_DTYPES = {
    "stable_id": "int64", "variant_id": "object", "chr": "object",
    "pos": "int64", "ref": "object", "alt": "object",
    "membership_count": "int64", "b_intercept": "float64",
    "a_diag": "float64",
    "n_nonmissing": "int64", "n_imputed": "int64",
}

COHORT_VARIANT_STATS_COLUMNS = [
    "stable_id", "b_intercept", "a_diag",
    "n_nonmissing", "n_imputed", "genotype_mean",
]

_COHORT_VARIANT_STATS_DTYPES = {
    "stable_id": "int64",
    "b_intercept": "float64",
    "a_diag": "float64",
    "n_nonmissing": "int64",
    "n_imputed": "int64",
    "genotype_mean": "float64",
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
        "cohort_stats_dir": d / "cohort_variant_stats",
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

    Returns (variants_df in file/row order, B (n_variants x 1), D (1 x 1)).
    Hard-fails if B's row count does not match the variant count. The
    intercept-only schema means B is a single column (per-variant dosage sum)
    and D is the 1x1 sample-count scalar; a legacy multi-column export is
    rejected rather than silently adjusted.
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
    if b_mat.shape[1] != 1 or d_mat.shape != (1, 1):
        raise ValueError(
            f"Expected B (n x 1) and D (1 x 1) for the intercept-only schema; "
            f"got B{b_mat.shape}, D{d_mat.shape}"
        )
    return variants, b_mat, d_mat


def merge_variant_table(
    table: pd.DataFrame,
    cohort_variants: pd.DataFrame,
    b_matrix: np.ndarray,
    next_stable_id: int,
) -> tuple[pd.DataFrame, dict[str, int], int]:
    """Add one cohort's per-variant data into the master table.

    Assigns a fresh stable_id to each unseen variant_id, sums the single B
    column (b_intercept = per-variant dosage sum), n_nonmissing, n_imputed,
    and bumps membership_count. a_diag is
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
                "membership_count": 1, "b_intercept": float(b[0]), "a_diag": 0.0,
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


def stage_merged_pairs(path: str | Path, incoming: np.ndarray, batch_rows: int) -> Path:
    """Sorted-merge-add ``incoming`` into a pair file, writing a staged temp.

    Streams the existing file in row batches (the side that grows with cohort
    count); the incoming array is assumed key-sorted and held in memory. The
    current file is left untouched; the merged result is written to a sibling
    ``.stage.<pid>`` file whose path is returned for a later atomic commit.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    incoming = incoming[np.argsort(_pair_key(incoming), kind="stable")]
    tmp = path.with_name(path.name + f".stage.{os.getpid()}")

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
    return tmp


def merge_pairs_into_file(path: str | Path, incoming: np.ndarray, batch_rows: int) -> None:
    """Stage then atomically commit a sorted-merge-add into a pair file."""
    path = Path(path)
    tmp = stage_merged_pairs(path, incoming, batch_rows)
    os.replace(tmp, path)


def stage_np_save(path: str | Path, arr: np.ndarray) -> tuple[Path, Path]:
    """Write ``arr`` to a staged ``.stage.<pid>`` temp; return (tmp, final)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".stage.{os.getpid()}.npy")
    np.save(tmp, arr)  # name already ends in .npy, so np.save writes exactly here
    return tmp, path


def stage_variant_table(precursor_dir: str | Path, df: pd.DataFrame) -> tuple[Path, Path]:
    """Write the variant table to a staged temp; return (tmp, final)."""
    path = precursor_paths(precursor_dir)["variants"]
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".stage.{os.getpid()}")
    df[VARIANT_COLUMNS].astype(_VARIANT_DTYPES).to_parquet(tmp, index=False)
    return tmp, path


def _safe_filename_component(value: str) -> str:
    """Return a conservative filename component derived from a cohort name."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._-")
    return safe or "cohort"


def cohort_variant_stats_path(
    precursor_dir: str | Path, cohort_index: int, study_name: str,
) -> Path:
    """Path for one cohort's retained per-variant contribution sidecar."""
    safe_name = _safe_filename_component(study_name)
    return (precursor_paths(precursor_dir)["cohort_stats_dir"] /
            f"{cohort_index:04d}_{safe_name}.parquet")


def stage_cohort_variant_stats(
    path: str | Path,
    cohort_variants: pd.DataFrame,
    b_matrix: np.ndarray,
    diag_acc: dict[int, float],
    id_map: dict[str, int],
) -> tuple[Path, Path]:
    """Stage one cohort's per-variant support, B, and A diagonal sidecar.

    The file deliberately stores stable IDs rather than repeating variant IDs;
    the master ``variants.parquet`` is the stable_id -> variant_id map. This is
    enough to reconstruct cohort-overlap centring later without inflating the
    precursor with repeated string keys.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f".stage.{os.getpid()}")

    stable_ids = np.array(
        [id_map[v] for v in cohort_variants["variant_id"]], dtype=np.int64,
    )
    a_diag = np.empty(len(stable_ids), dtype=np.float64)
    missing: list[str] = []
    for idx, sid in enumerate(stable_ids):
        val = diag_acc.get(int(sid))
        if val is None:
            missing.append(str(cohort_variants.iloc[idx]["variant_id"]))
            a_diag[idx] = np.nan
        else:
            a_diag[idx] = val
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(
            f"Missing A diagonal entries for {len(missing)} variants while "
            f"writing cohort variant stats; first missing: {preview}"
        )

    df = pd.DataFrame({
        "stable_id": stable_ids,
        "b_intercept": b_matrix[:, 0].astype(np.float64),
        "a_diag": a_diag,
        "n_nonmissing": cohort_variants["n_nonmissing"].to_numpy(np.int64),
        "n_imputed": cohort_variants["n_imputed"].to_numpy(np.int64),
        "genotype_mean": cohort_variants["genotype_mean"].to_numpy(np.float64),
    })
    df[COHORT_VARIANT_STATS_COLUMNS].astype(
        _COHORT_VARIANT_STATS_DTYPES,
    ).to_parquet(tmp, index=False)
    return tmp, path


def _commit_staged(staged: list[tuple[Path, Path]]) -> None:
    """Atomically rename each staged temp into its final path (commit burst)."""
    for tmp, final in staged:
        os.replace(tmp, final)


def _clean_stale_temps(precursor_dir: Path) -> None:
    """Remove orphaned staging temps left by a crashed prior run (lock-guarded)."""
    paths = precursor_paths(precursor_dir)
    for base in (paths["root"], paths["pairs_dir"], paths["cohort_stats_dir"]):
        if base.is_dir():
            for pattern in ("*.stage.*", "*.tmp.*"):
                for stale in base.glob(pattern):
                    stale.unlink(missing_ok=True)


def _incomplete_cohorts(pm: dict) -> list[str]:
    """Study names of any cohort entry not in the committed state."""
    return [c["study_name"] for c in pm.get("cohorts", [])
            if c.get("status") != "committed"]


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
    verify_checksums: bool = True,
) -> None:
    """Add one cohort's A/B/D into the precursor (see spec accumulate phase).

    Integrity gate: unless ``verify_checksums`` is disabled, the cohort's
    ``checksums.json`` is re-verified first, so a corrupted/truncated upload
    fails loudly before any precursor state is touched.

    Crash safety: the per-cohort mutations are *staged* to sibling temp files
    while the live precursor is untouched, then a write-ahead ``in_progress``
    manifest entry is recorded, then all temps are renamed in (commit burst),
    then the entry is flipped to ``committed``. A crash during the commit window
    leaves the entry ``in_progress``; a later ``accumulate``/``finalise`` then
    refuses to proceed (rather than silently double-counting) until the operator
    restores or rebuilds. ``force=True`` overrides both the duplicate-study and
    incomplete-entry guards.
    """
    cohort_dir = Path(cohort_dir)
    precursor_dir = Path(precursor_dir)
    if verify_checksums:
        ld_checksums.verify_cohort_checksums(cohort_dir)
    cohort_manifest = json.loads((cohort_dir / "manifest.json").read_text())
    study_name = cohort_manifest["study_name"]
    contract = extract_contract(cohort_manifest)

    pm = read_precursor_manifest(precursor_dir)
    if pm is None:
        pm = {
            "schema_version": PRECURSOR_SCHEMA_VERSION,
            "cohort_variant_stats_schema_version": COHORT_VARIANT_STATS_SCHEMA_VERSION,
            "panel_specification_version": DEFAULT_PANEL_SPEC_VERSION,
            "contract": contract, "n_cohorts": 0, "cohorts": [],
            "next_stable_id": 0,
        }
        precursor_dir.mkdir(parents=True, exist_ok=True)
    else:
        validate_contract(pm["contract"], contract)
        pm.setdefault(
            "cohort_variant_stats_schema_version",
            COHORT_VARIANT_STATS_SCHEMA_VERSION,
        )
        incomplete = _incomplete_cohorts(pm)
        if incomplete and not force:
            raise RuntimeError(
                f"Precursor has incomplete (in-progress) accumulate(s) for "
                f"{incomplete}; it may be partially updated. Restore from backup "
                "or rebuild, then retry (or pass force=True to override)."
            )
        if any(c["study_name"] == study_name for c in pm["cohorts"]) and not force:
            raise ValueError(
                f"Study '{study_name}' already accumulated into this precursor; "
                "rebuild or pass force=True to override."
            )

    lock = _acquire_lock(precursor_dir)
    try:
        _clean_stale_temps(precursor_dir)
        variants, b_mat, d_mat = read_cohort_assets(cohort_dir)
        table = read_variant_table(precursor_dir)
        cohort_index = int(pm["n_cohorts"])
        table, id_map, next_id = merge_variant_table(
            table, variants, b_mat, pm["next_stable_id"])

        # cohort per-chrom canonical arrays for index -> (pos, sid) mapping
        cohort_chrom = {}
        for chrom, sub in variants.groupby("chr", sort=False):
            cohort_chrom[chrom] = {
                "pos": sub["pos"].to_numpy(np.int64),
                "sid": np.array([id_map[v] for v in sub["variant_id"]], dtype=np.int64),
            }

        # --- compute + stage everything; live precursor stays untouched ---
        staged: list[tuple[Path, Path]] = []
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
                pair_final = precursor_paths(precursor_dir)["pairs_dir"] / f"chr{chrom}.parquet"
                pair_tmp = stage_merged_pairs(pair_final, incoming, batch_rows=pair_batch_rows)
                staged.append((pair_tmp, pair_final))

        # apply diagonal accumulation (in memory) and stage the variant table
        if diag_acc:
            add = table["stable_id"].map(lambda s: diag_acc.get(int(s), 0.0))
            table["a_diag"] = table["a_diag"] + add
        stats_final = cohort_variant_stats_path(
            precursor_dir, cohort_index, study_name,
        )
        staged.append(stage_cohort_variant_stats(
            stats_final, variants, b_mat, diag_acc, id_map,
        ))
        staged.append(stage_variant_table(precursor_dir, table))

        d_path = precursor_paths(precursor_dir)["d"]
        existing_d = np.load(d_path) if d_path.is_file() else np.zeros((1, 1))
        staged.append(stage_np_save(d_path, existing_d + d_mat))

        # --- write-ahead intent, commit, then mark committed (tight window) ---
        pm["next_stable_id"] = next_id
        pm["n_cohorts"] += 1
        pm["cohorts"].append({
            "study_name": study_name,
            "accumulated_at_utc": datetime.now(timezone.utc).isoformat(),
            "n_variants": int(len(variants)),
            "d_intercept": float(d_mat[0, 0]),
            "cohort_variant_stats": stats_final.relative_to(precursor_dir).as_posix(),
            "status": "in_progress",
        })
        write_precursor_manifest(precursor_dir, pm)

        _commit_staged(staged)

        pm["cohorts"][-1]["status"] = "committed"
        write_precursor_manifest(precursor_dir, pm)
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
    b = kept[["b_intercept"]].to_numpy()
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


class _PairCursor:
    """Sequential reader over a position-sorted A_pairs/chr*.parquet file.

    On construction, each stored pair's endpoints are mapped sid -> pooled index;
    pairs with an endpoint outside the intersection (sid not in the map) are
    dropped. Each surviving pair is normalised to the upper triangle by pooled
    index (lo <= hi), so co-located variants (same position, differing ref/alt)
    whose sid order disagrees with the pooled/canonical order still scatter into
    the correct upper-triangular cell. Pairs are sorted by the pooled row index
    `lo` so fill_block can stream them in chunk order with a safe early break.
    Indices are made chromosome-local by subtracting `chrom_first_pooled` so they
    match finalise's per-chromosome row_start/row_stop (which reset to 0).
    """
    def __init__(self, path, sid_to_pooled, chrom_first_pooled=0):
        if chrom_first_pooled is None:
            chrom_first_pooled = 0
        if Path(path).is_file():
            df = pq.read_table(path).to_pandas()
            pi = df["sid_i"].map(sid_to_pooled)
            pj = df["sid_j"].map(sid_to_pooled)
            keep = (pi.notna() & pj.notna()).to_numpy()
            pi = pi.to_numpy()[keep].astype(np.int64)
            pj = pj.to_numpy()[keep].astype(np.int64)
            val = df["value"].to_numpy()[keep]
            # chromosome-local pooled indices (finalise's row_start/row_stop reset
            # to 0 per chromosome, but pooled_index is global)
            lo = np.minimum(pi, pj) - chrom_first_pooled
            hi = np.maximum(pi, pj) - chrom_first_pooled
            order = np.argsort(lo, kind="stable")   # monotonic row index for streaming
            self._lo = lo[order]
            self._hi = hi[order]
            self._val = val[order]
            self._n = int(self._lo.size)
        else:
            self._lo = self._hi = self._val = None
            self._n = 0
        self._cursor = 0

    def fill_block(self, block, row_start, row_stop, col_stop, idx):
        # col_stop/idx kept for signature compatibility with the caller; the
        # window bound is guaranteed by construction (hi < local_stops[lo]).
        while self._cursor < self._n:
            r = int(self._lo[self._cursor])
            if r >= row_stop:
                break
            if r < row_start:
                self._cursor += 1
                continue
            c = int(self._hi[self._cursor])
            block[r - row_start, c - row_start] += self._val[self._cursor]
            self._cursor += 1


def _reserve_panel_version_dir(panel_root: Path) -> tuple[Path, Path]:
    """Return (work_dir, final_dir) for the next immutable panel_vN output."""
    panel_root.mkdir(parents=True, exist_ok=True)
    used = []
    for child in panel_root.iterdir():
        match = re.match(r"^panel_v([0-9]+)$", child.name)
        if match and child.is_dir():
            used.append(int(match.group(1)))
    version = max(used, default=0) + 1
    final_dir = panel_root / f"panel_v{version}"
    work_dir = panel_root / f".{final_dir.name}.tmp.{os.getpid()}"
    work_dir.mkdir(parents=True, exist_ok=False)
    return work_dir, final_dir


def finalise(precursor_dir, panel_dir, r_writer, maf_threshold, min_adj_diag,
             block_size, max_dense_gb=1.0, min_cohorts=None):
    """Resolve intersection, adjust, convert to R, write a versioned panel."""
    precursor_dir, panel_dir = Path(precursor_dir), Path(panel_dir)
    pm = read_precursor_manifest(precursor_dir)
    if pm is None or pm["n_cohorts"] == 0:
        raise ValueError("Cannot finalise an empty precursor")
    incomplete = _incomplete_cohorts(pm)
    if incomplete:
        raise RuntimeError(
            f"Precursor has incomplete (in-progress) accumulate(s) for "
            f"{incomplete}; refusing to finalise a possibly partial precursor. "
            "Restore from backup or rebuild first."
        )
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
    b_all = surv[["b_intercept"]].to_numpy()
    w_all = b_all @ d_inv
    adj_diag_all = surv["a_adj_diag"].to_numpy()
    sid_to_pooled = dict(zip(surv["stable_id"].to_numpy(), surv["pooled_index"].to_numpy()))
    pos_all = surv["pos"].to_numpy(np.int64)

    panel_work_dir, panel_final_dir = _reserve_panel_version_dir(panel_dir)
    r_root = panel_work_dir / "R_blocks"

    for chrom, sub in surv.groupby("chr", sort=False):
        idx = sub["pooled_index"].to_numpy(np.int64)
        local_pos = pos_all[idx]
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

    _write_panel_artefacts(panel_work_dir, pm, surv, dropped, d_matrix, d_rank,
                           maf_threshold, min_adj_diag, panel_final_dir.name)
    panel_work_dir.rename(panel_final_dir)
    return panel_final_dir


def _write_panel_artefacts(panel_dir, pm, surv, dropped, d_matrix, d_rank,
                           maf_threshold, min_adj_diag, panel_version):
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
        "module": "15_aggregate", "panel_version": panel_version,
        "panel_specification_version": pm["panel_specification_version"],
        "covariate_schema": pm["contract"]["schema_id"], "n_cohorts": pm["n_cohorts"],
        "cohorts": [c["study_name"] for c in pm["cohorts"]],
        "n_variants_panel": int(len(surv)), "n_variants_dropped": int(len(dropped)),
        "d_rank": d_rank, "d_condition_number": float(np.linalg.cond(d_matrix)),
        "maf_threshold": maf_threshold, "min_adj_diag": min_adj_diag,
        "regularisation": "none", "r_dtype": "float64",
    }
    (panel_dir / "pooled_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    (panel_dir / "qc_report.txt").write_text(
        f"Section 15 aggregate panel built: {len(surv)} variants, {len(dropped)} dropped, "
        f"{pm['n_cohorts']} cohorts, D rank {d_rank}.\n")
