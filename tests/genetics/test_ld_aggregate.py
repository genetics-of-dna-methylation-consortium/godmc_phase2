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
