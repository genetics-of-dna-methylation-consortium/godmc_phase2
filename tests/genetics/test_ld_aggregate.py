# tests/genetics/test_ld_aggregate.py
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import ld_aggregate as agg
import ld_checksums as ck


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
        [[0, "1:100:G:A", "1", 100, "G", "A", 1, 3.0, 5.0, 4, 0]],
        columns=agg.VARIANT_COLUMNS,
    )
    agg.write_variant_table(tmp_path, df)
    back = agg.read_variant_table(tmp_path)
    pd.testing.assert_frame_equal(back, df)


def _cohort_manifest():
    return {
        "study_name": "cohortA",
        "genome_build": "GRCh37",
        "covariate_schema": {
            "schema_id": "intercept_only",
            "matrix_columns": ["intercept"],
            "sex_factor_recode": {},
        },
        "variant_index": {"schema_version": "v0.3-with-genotype-stats"},
        "A_blocks": {"radius_bp": 1_000_000, "block_size": 4096},
    }


def test_extract_contract_pulls_expected_fields():
    c = agg.extract_contract(_cohort_manifest())
    assert c == {
        "genome_build": "GRCh37",
        "schema_id": "intercept_only",
        "matrix_columns": ["intercept"],
        "sex_factor_recode": {},
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
    b = [[4.0], [2.0]]
    d = [[7.0]]
    out = _write_cohort_dir(tmp_path, rows, b, d)

    variants, b_mat, d_mat = agg.read_cohort_assets(out)
    assert list(variants["variant_id"]) == ["1:100:G:A", "1:200:C:T"]
    assert variants["pos"].tolist() == [100, 200]
    np.testing.assert_array_equal(b_mat, np.asarray(b))
    np.testing.assert_array_equal(d_mat, np.asarray(d))
    assert b_mat.shape == (2, 1)


def test_read_cohort_assets_rejects_b_row_mismatch(tmp_path):
    rows = [("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0)]
    out = _write_cohort_dir(tmp_path, rows, [[1.0], [4.0]], [[1.0]])
    with pytest.raises(ValueError, match="rows"):
        agg.read_cohort_assets(out)


def test_read_cohort_assets_rejects_non_intercept_only_shapes(tmp_path):
    # The intercept-only schema requires B (n x 1) and D (1 x 1); a legacy
    # 3-column export (intercept+age+sex) must hard-fail rather than be adjusted.
    rows = [("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0)]
    out = _write_cohort_dir(tmp_path, rows, [[4.0, 90.0, 1.5]], np.eye(3))
    with pytest.raises(ValueError, match="1"):
        agg.read_cohort_assets(out)


def test_variant_columns_have_no_age_sex():
    assert "b_age" not in agg.VARIANT_COLUMNS
    assert "b_sex" not in agg.VARIANT_COLUMNS
    assert "b_intercept" in agg.VARIANT_COLUMNS


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
    b = np.array([[4.0], [2.0]])

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
    b2 = np.array([[6.0], [3.0]])
    table, id_map2, next_id = agg.merge_variant_table(table, cohort2, b2, next_stable_id=next_id)

    assert next_id == 3
    assert id_map2 == {"1:100:G:A": 0, "1:300:A:G": 2}
    shared = table.set_index("variant_id").loc["1:100:G:A"]
    assert shared["membership_count"] == 2
    assert shared["b_intercept"] == 10.0       # 4 + 6
    assert shared["n_nonmissing"] == 9         # 4 + 5
    assert table.set_index("variant_id").loc["1:300:A:G"]["membership_count"] == 1


def test_extract_chunk_entries_splits_diag_and_offdiag():
    # 3 variants on one chrom, positions 100/200/2_000_000, radius 1Mb.
    positions = np.array([100, 200, 2_000_000], dtype=np.int64)
    sids = np.array([10, 11, 12], dtype=np.int64)
    radius = 1_000_000
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


def _build_synthetic_cohort(tmp_path, study_name="cohortA"):
    """Build a tiny 2-variant chr1 cohort dir (variants/B/D/manifest) for accumulate tests."""
    cohort = tmp_path / study_name
    cohort.mkdir()
    rows = [("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0),
            ("1", 200, "C", "T", "1:200:C:T", 4, 0, 0.5)]
    with gzip.open(cohort / "variants.tsv.gz", "wt") as fh:
        fh.write("chr\tpos\tref\talt\tvariant_id\tn_nonmissing\tn_imputed\tgenotype_mean\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    np.save(cohort / "B.npy", np.array([[4.0], [2.0]]))
    np.save(cohort / "D.npy", np.array([[4.0]]))
    manifest = {
        "study_name": study_name, "genome_build": "GRCh37",
        "covariate_schema": {"schema_id": "intercept_only",
            "matrix_columns": ["intercept"], "sex_factor_recode": {}},
        "variant_index": {"schema_version": "v0.3-with-genotype-stats"},
        "A_blocks": {"radius_bp": 1_000_000, "block_size": 8,
            "chromosomes": {"1": {"chunks": [
                {"row_start": 0, "row_stop": 2, "column_start": 0,
                 "column_stop": 2, "directory": "A_blocks/chr1/chunk_000000"}]}}},
    }
    (cohort / "manifest.json").write_text(json.dumps(manifest))
    ck.write_cohort_checksums(cohort)
    return cohort


# symmetric dense A for the single chunk: X X^T (diag 20,10; offdiag 8)
def _fake_reader(chunk_dir):
    return np.array([[20.0, 8.0], [8.0, 10.0]])


def test_accumulate_one_cohort_builds_precursor(tmp_path):
    cohort = _build_synthetic_cohort(tmp_path)
    precursor = tmp_path / "precursor"
    agg.accumulate(cohort, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)

    pm = agg.read_precursor_manifest(precursor)
    assert pm["n_cohorts"] == 1
    assert pm["cohorts"][0]["study_name"] == "cohortA"
    table = agg.read_variant_table(precursor).set_index("variant_id")
    assert table.loc["1:100:G:A", "a_diag"] == 20.0
    assert table.loc["1:200:C:T", "a_diag"] == 10.0
    pairs = pq.read_table(precursor / "A_pairs" / "chr1.parquet").to_pandas()
    assert pairs["value"].tolist() == [8.0]
    np.testing.assert_array_equal(np.load(precursor / "D.npy"),
                                  np.load(cohort / "D.npy"))


def test_accumulate_rejects_duplicate_study(tmp_path):
    cohort = _build_synthetic_cohort(tmp_path)
    precursor = tmp_path / "precursor"
    agg.accumulate(cohort, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)
    with pytest.raises(ValueError, match="already accumulated"):
        agg.accumulate(cohort, precursor, chunk_reader=_fake_reader)


def test_accumulate_verifies_checksums_and_rejects_tampered_cohort(tmp_path):
    cohort = _build_synthetic_cohort(tmp_path)
    # tamper with B.npy AFTER checksums were written by the builder
    np.save(cohort / "B.npy", np.array([[99.0], [99.0]]))
    precursor = tmp_path / "precursor"
    with pytest.raises(ck.ChecksumError, match="B.npy"):
        agg.accumulate(cohort, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)
    # nothing was committed
    assert agg.read_precursor_manifest(precursor) is None


def test_accumulate_rejects_cohort_without_checksums(tmp_path):
    cohort = _build_synthetic_cohort(tmp_path)
    (cohort / ck.CHECKSUM_FILENAME).unlink()
    precursor = tmp_path / "precursor"
    with pytest.raises(ck.ChecksumError, match=ck.CHECKSUM_FILENAME):
        agg.accumulate(cohort, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)


def test_accumulate_marks_cohort_committed(tmp_path):
    cohort = _build_synthetic_cohort(tmp_path)
    precursor = tmp_path / "precursor"
    agg.accumulate(cohort, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)
    pm = agg.read_precursor_manifest(precursor)
    assert pm["cohorts"][0]["status"] == "committed"


def test_accumulate_refuses_when_in_progress_entry_present(tmp_path):
    cohort_a = _build_synthetic_cohort(tmp_path, study_name="cohortA")
    precursor = tmp_path / "precursor"
    agg.accumulate(cohort_a, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)
    # simulate a prior crashed accumulate: an entry stuck in_progress
    pm = agg.read_precursor_manifest(precursor)
    pm["cohorts"][0]["status"] = "in_progress"
    agg.write_precursor_manifest(precursor, pm)

    cohort_b = _build_synthetic_cohort(tmp_path, study_name="cohortB")
    with pytest.raises(RuntimeError, match="incomplete"):
        agg.accumulate(cohort_b, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)


def test_accumulate_crash_at_commit_leaves_in_progress_and_no_partial_data(tmp_path, monkeypatch):
    cohort_a = _build_synthetic_cohort(tmp_path, study_name="cohortA")
    precursor = tmp_path / "precursor"
    agg.accumulate(cohort_a, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)
    before = agg.read_variant_table(precursor).set_index("variant_id")["b_intercept"].to_dict()
    d_before = np.load(precursor / "D.npy").copy()

    cohort_b = _build_synthetic_cohort(tmp_path, study_name="cohortB")

    def boom(staged):
        raise RuntimeError("disk died mid-commit")
    monkeypatch.setattr(agg, "_commit_staged", boom)
    with pytest.raises(RuntimeError, match="mid-commit"):
        agg.accumulate(cohort_b, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)

    # write-ahead intent recorded but data NOT applied (no double-count)
    pm = agg.read_precursor_manifest(precursor)
    assert pm["cohorts"][-1]["study_name"] == "cohortB"
    assert pm["cohorts"][-1]["status"] == "in_progress"
    after = agg.read_variant_table(precursor).set_index("variant_id")["b_intercept"].to_dict()
    assert after == before                                    # data untouched
    np.testing.assert_array_equal(np.load(precursor / "D.npy"), d_before)

    # and a subsequent normal accumulate refuses until the operator intervenes
    monkeypatch.undo()
    cohort_c = _build_synthetic_cohort(tmp_path, study_name="cohortC")
    with pytest.raises(RuntimeError, match="incomplete"):
        agg.accumulate(cohort_c, precursor, chunk_reader=_fake_reader, pair_batch_rows=16)


def test_finalise_refuses_in_progress_precursor(tmp_path):
    precursor = _build_two_cohort_precursor(tmp_path)
    pm = agg.read_precursor_manifest(precursor)
    pm["cohorts"][0]["status"] = "in_progress"
    agg.write_precursor_manifest(precursor, pm)
    with pytest.raises(RuntimeError, match="incomplete"):
        agg.finalise(precursor, tmp_path / "panel", r_writer=lambda *a, **k: None,
                     maf_threshold=0.0, min_adj_diag=-1e9, block_size=8)


def test_resolve_intersection_keeps_full_membership():
    table = pd.DataFrame({
        "stable_id": [0, 1, 2], "variant_id": ["1:200:C:T", "1:100:G:A", "2:50:A:C"],
        "chr": ["1", "1", "2"], "pos": [200, 100, 50], "ref": ["C", "G", "A"],
        "alt": ["T", "A", "C"], "membership_count": [2, 1, 2],
        "b_intercept": [0.0]*3,
        "a_diag": [0.0]*3, "n_nonmissing": [0]*3, "n_imputed": [0]*3,
    })
    kept = agg.resolve_intersection(table, n_cohorts=2, min_cohorts=2)
    # variant 1 (membership 1) dropped; survivors sorted by chr,pos
    assert kept["variant_id"].tolist() == ["1:200:C:T", "2:50:A:C"]
    assert kept["pooled_index"].tolist() == [0, 1]


def test_apply_filters_drops_rare_and_unstable():
    # D[0,0] = N = 100. b_intercept 30 -> f=0.15 keep; b_intercept 1 -> f=0.005 drop
    d = np.array([[100.0]])
    kept = pd.DataFrame({
        "stable_id": [0, 1], "variant_id": ["1:100:G:A", "1:200:C:T"],
        "chr": ["1", "1"], "pos": [100, 200], "ref": ["G", "C"], "alt": ["A", "T"],
        "b_intercept": [30.0, 1.0],
        "a_diag": [50.0, 50.0], "pooled_index": [0, 1],
    })
    surv, dropped = agg.apply_filters(kept, d, maf_threshold=0.01, min_adj_diag=0.0)
    assert surv["variant_id"].tolist() == ["1:100:G:A"]
    assert surv["pooled_index"].tolist() == [0]                 # re-indexed
    assert "1:200:C:T" in dropped["variant_id"].tolist()
    assert dropped.set_index("variant_id").loc["1:200:C:T", "reason"] == "maf_below_threshold"
    assert "a_adj_diag" in surv.columns


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


def _build_two_cohort_precursor(tmp_path):
    """Build a 2-cohort precursor by accumulating two synthetic cohorts (overlapping
    3 variants on chr1) with a fake chunk_reader returning each cohort's X X^T.
    Returns the precursor path. Also used to anchor finalise tests."""
    variant_meta = [("1", 100, "G", "A", "1:100:G:A"),
                    ("1", 200, "C", "T", "1:200:C:T"),
                    ("1", 300, "A", "G", "1:300:A:G")]
    precursor = tmp_path / "precursor"

    def build_and_accumulate(name, n_samples, seed):
        r = np.random.default_rng(seed)
        X = r.integers(0, 3, size=(n_samples, 3)).astype(float)
        C = np.ones((n_samples, 1))
        A = X.T @ X
        B = X.T @ C
        D = C.T @ C
        cohort = tmp_path / name
        cohort.mkdir()
        with gzip.open(cohort / "variants.tsv.gz", "wt") as fh:
            fh.write("chr\tpos\tref\talt\tvariant_id\tn_nonmissing\tn_imputed\tgenotype_mean\n")
            for col, (c, p, ref, alt, vid) in enumerate(variant_meta):
                fh.write(f"{c}\t{p}\t{ref}\t{alt}\t{vid}\t{n_samples}\t0\t{X[:, col].mean()}\n")
        np.save(cohort / "B.npy", B)
        np.save(cohort / "D.npy", D)
        manifest = {
            "study_name": name, "genome_build": "GRCh37",
            "covariate_schema": {"schema_id": "intercept_only",
                "matrix_columns": ["intercept"], "sex_factor_recode": {}},
            "variant_index": {"schema_version": "v0.3-with-genotype-stats"},
            "A_blocks": {"radius_bp": 1_000_000, "block_size": 8,
                "chromosomes": {"1": {"chunks": [
                    {"row_start": 0, "row_stop": 3, "column_start": 0,
                     "column_stop": 3, "directory": "A_blocks/chr1/chunk_000000"}]}}},
        }
        (cohort / "manifest.json").write_text(json.dumps(manifest))
        ck.write_cohort_checksums(cohort)
        agg.accumulate(cohort, precursor, chunk_reader=lambda cd: A, pair_batch_rows=64)

    build_and_accumulate("cohort_a", 30, 1)
    build_and_accumulate("cohort_b", 40, 2)
    return precursor


def test_finalise_writes_panel_and_calls_writer(tmp_path):
    precursor = _build_two_cohort_precursor(tmp_path)
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


def test_paircursor_upper_triangle_and_no_drop_across_chunks(tmp_path):
    # sid -> pooled with co-located disorder: sid 0 is canonically LATER (pooled 1),
    # sid 1 is canonically EARLIER (pooled 0); plus sid 2 -> pooled 2, sid 3 -> pooled 3.
    sid_to_pooled = {0: 1, 1: 0, 2: 2, 3: 3}
    # stored pairs (pos_i, sid_i, pos_j, sid_j, value), key-sorted by (pos_i, sid_i, ...):
    #  (100, 0, 200, 2, 5.0) -> pooled (1,2) upper
    #  (100, 1, 200, 3, 7.0) -> pooled (0,3) upper
    #  (100, 0, 100, 1, 9.0) -> co-located, pooled (1,0) -> must normalise to (0,1)
    rows = [(100, 0, 200, 2, 5.0), (100, 1, 200, 3, 7.0), (100, 0, 100, 1, 9.0)]
    arr = np.empty(len(rows), dtype=agg.PAIR_DTYPE)
    for i, r in enumerate(rows):
        arr[i] = r
    path = tmp_path / "chr1.parquet"
    agg.merge_pairs_into_file(path, arr, batch_rows=8)

    cur = agg._PairCursor(path, sid_to_pooled, chrom_first_pooled=0)
    # 4 pooled rows (0..3); simulate TWO chunks split at row 2
    n = 4
    full = np.zeros((n, n))
    # chunk 0: rows [0,2), cols up to 4
    b0 = np.zeros((2, n)); cur.fill_block(b0, 0, 2, n, np.arange(n))
    full[0:2, :] += b0
    # chunk 1: rows [2,4)
    b1 = np.zeros((2, n)); cur.fill_block(b1, 2, 4, n, np.arange(n))
    full[2:4, :] += b1

    # expected upper-triangle placements:
    #  (1,2)=5.0, (0,3)=7.0, (0,1)=9.0  (the co-located pair normalised to upper)
    assert full[1, 2] == 5.0
    assert full[0, 3] == 7.0
    assert full[0, 1] == 9.0            # was lower-triangle (1,0) before the fix
    # nothing in the lower triangle
    assert np.allclose(np.tril(full, -1), 0.0)
    # total mass conserved (no dropped pair across the chunk boundary)
    assert full.sum() == 21.0


def test_finalise_multi_chromosome_keeps_offdiagonals(tmp_path):
    # Two cohorts, 2 variants on chr1 + 2 on chr2 (full overlap). Catches the
    # global-vs-chromosome-local pooled-index bug: chr2 off-diagonal pairs must
    # NOT be dropped, and chr2 R must match the direct per-chromosome reference.
    chroms = ["1", "2"]
    variant_meta = [
        ("1", 100, "G", "A", "1:100:G:A"),
        ("1", 200, "C", "T", "1:200:C:T"),
        ("2", 100, "A", "C", "2:100:A:C"),
        ("2", 200, "T", "G", "2:200:T:G"),
    ]
    chrom_cols = {"1": [0, 1], "2": [2, 3]}
    precursor = tmp_path / "precursor"
    cohort_X, cohort_C = {}, {}

    def build_and_accumulate(name, n_samples, seed):
        r = np.random.default_rng(seed)
        X = r.integers(0, 3, size=(n_samples, 4)).astype(float)
        C = np.ones((n_samples, 1))
        cohort_X[name], cohort_C[name] = X, C
        cohort = tmp_path / name
        cohort.mkdir()
        with gzip.open(cohort / "variants.tsv.gz", "wt") as fh:
            fh.write("chr\tpos\tref\talt\tvariant_id\tn_nonmissing\tn_imputed\tgenotype_mean\n")
            for col, (c, p, ref, alt, vid) in enumerate(variant_meta):
                fh.write(f"{c}\t{p}\t{ref}\t{alt}\t{vid}\t{n_samples}\t0\t{X[:, col].mean()}\n")
        np.save(cohort / "B.npy", X.T @ C)
        np.save(cohort / "D.npy", C.T @ C)
        manifest = {
            "study_name": name, "genome_build": "GRCh37",
            "covariate_schema": {"schema_id": "intercept_only",
                "matrix_columns": ["intercept"], "sex_factor_recode": {}},
            "variant_index": {"schema_version": "v0.3-with-genotype-stats"},
            "A_blocks": {"radius_bp": 1_000_000, "block_size": 8, "chromosomes": {
                c: {"chunks": [{"row_start": 0, "row_stop": 2, "column_start": 0,
                    "column_stop": 2, "directory": f"A_blocks/chr{c}/chunk_000000"}]}
                for c in chroms}},
        }
        (cohort / "manifest.json").write_text(json.dumps(manifest))
        ck.write_cohort_checksums(cohort)

        def reader(chunk_dir, _X=X):
            cols = chrom_cols["1"] if "chr1" in str(chunk_dir) else chrom_cols["2"]
            sub = _X[:, cols]
            return sub.T @ sub

        agg.accumulate(cohort, precursor, chunk_reader=reader, pair_batch_rows=64)

    build_and_accumulate("cohort_a", 30, 1)
    build_and_accumulate("cohort_b", 40, 2)

    panel = tmp_path / "panel"
    captured = {}
    def capture_writer(dense, starts, stops, out_dir, block_size):
        chrom = "1" if "chr1" in str(out_dir) else "2"
        captured.setdefault(chrom, []).append(np.asarray(dense).copy())
        Path(out_dir).mkdir(parents=True, exist_ok=True)

    agg.finalise(precursor, panel, r_writer=capture_writer,
                 maf_threshold=0.0, min_adj_diag=-1e9, block_size=8, max_dense_gb=1.0)

    Xf = np.vstack([cohort_X["cohort_a"], cohort_X["cohort_b"]])
    Cf = np.vstack([cohort_C["cohort_a"], cohort_C["cohort_b"]])
    B = Xf.T @ Cf
    d_inv = np.linalg.inv(Cf.T @ Cf)
    for chrom, cols in chrom_cols.items():
        Xc = Xf[:, cols]
        Bc = B[cols, :]
        A_adj = Xc.T @ Xc - Bc @ d_inv @ Bc.T
        diag = np.diag(A_adj)
        R_ref = A_adj / np.sqrt(np.outer(diag, diag))
        assert chrom in captured, f"no R block emitted for chr{chrom}"
        R_got = captured[chrom][0]
        np.testing.assert_allclose(np.triu(R_got), np.triu(R_ref), rtol=1e-9, atol=1e-12)
        assert abs(np.triu(R_got, 1).sum()) > 0    # off-diagonal present (bug zeroed chr2's)
