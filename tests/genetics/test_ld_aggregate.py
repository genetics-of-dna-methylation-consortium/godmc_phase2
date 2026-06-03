# tests/genetics/test_ld_aggregate.py
import gzip
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
