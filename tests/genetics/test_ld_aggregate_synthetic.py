"""Synthetic two-cohort federated-equivalence harness for section 15.

Runs ``ld_prepare_stats.py`` against two independent synthetic cohorts and
checks that summing the cohort sufficient statistics matches a direct
stacked-residualisation reference computed in NumPy.

This module exercises the full section-15b pipeline end to end: real 15a
(``ld_prepare_stats.py``) -> ``accumulate`` -> ``finalise``, then reads the
panel ``R_blocks`` back and asserts equality against a direct
stacked-residualisation NumPy reference. It covers ``D``/``B`` federated
equivalence, ``A_adj``/``R`` equivalence (masked to the LD radius window),
accumulate order-independence, variant-set intersection, and the
duplicate-study accumulation guard.
"""

import gzip
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from synthetic_plink import write_synthetic_plink

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "resources" / "genetics" / "ld_prepare_stats.py"
AGG = REPO_ROOT / "resources" / "genetics" / "ld_aggregate_stats.py"

pytest.importorskip("hail", reason="hail not installed")


def _accumulate(cohort_out, precursor, tmp_path):
    r = subprocess.run([sys.executable, str(AGG), "--mode", "accumulate",
        "--cohort-dir", str(cohort_out), "--precursor-dir", str(precursor),
        "--log-file", str(tmp_path / "agg.log")], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def _finalise(precursor, panel, tmp_path):
    r = subprocess.run([sys.executable, str(AGG), "--mode", "finalise",
        "--precursor-dir", str(precursor), "--panel-dir", str(panel),
        "--maf-threshold", "0.0", "--min-adj-diag=-1e9",
        "--a-block-size", "8", "--log-file", str(tmp_path / "fin.log")],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def _read_r_full(panel, chrom, n):
    """Reassemble the full pooled R matrix for a chromosome from its chunks."""
    from hail.linalg import BlockMatrix
    out = np.zeros((n, n))
    chunk_root = panel / "R_blocks" / f"chr{chrom}"
    for c in sorted(chunk_root.glob("chunk_*")):
        block = BlockMatrix.read(str(c)).to_numpy()
        out[0:block.shape[0], 0:block.shape[1]] += block   # single-chunk chr22 here
    return out


SHARED_VARIANTS = [
    ("22", 1_000_000, "A", "G"),
    ("22", 2_000_000, "C", "T"),
    ("22", 3_000_000, "T", "A"),
]

COHORT_A_CALLS = [
    [0, 1, 2, 0],
    [1, 1, 0, 2],
    [2, 0, 1, 1],
]
COHORT_A_SAMPLES = ["A_IID1", "A_IID2", "A_IID3", "A_IID4"]
COHORT_A_COVARIATES = [
    ("A_IID1", 30, "M"),
    ("A_IID2", 45, "F"),
    ("A_IID3", 60, "M"),
    ("A_IID4", 55, "F"),
]

COHORT_B_CALLS = [
    [1, 0, 2, 1, 0],
    [2, 1, 0, 1, 2],
    [0, 2, 1, 0, 1],
]
COHORT_B_SAMPLES = ["B_IID1", "B_IID2", "B_IID3", "B_IID4", "B_IID5"]
COHORT_B_COVARIATES = [
    ("B_IID1", 25, "M"),
    ("B_IID2", 40, "F"),
    ("B_IID3", 50, "M"),
    ("B_IID4", 65, "F"),
    ("B_IID5", 35, "M"),
]


def _covariate_matrix(rows: list[tuple[str, int, str]]) -> np.ndarray:
    # Section 15 is intercept-only (grand-mean centring): Age/Sex in the cohort
    # covariate file are deliberately ignored, so the reference C is a column of
    # ones. With this C, A_adj is the mean-centred cross-product and R reduces to
    # the plain Pearson correlation of stacked dosages -> the GoDMC mQTL estimand.
    return np.ones((len(rows), 1), dtype=np.float64)


def _write_covariates_file(path: Path, rows: list[tuple[str, int, str]]) -> None:
    body = "IID\tAge_numeric\tSex_factor\n" + "\n".join(
        f"{iid}\t{age}\t{sex}" for iid, age, sex in rows
    )
    path.write_text(body + "\n")


def _read_variant_ids(out_dir: Path) -> list[str]:
    with gzip.open(out_dir / "variants.tsv.gz", "rt") as handle:
        lines = handle.read().splitlines()
    return [line.split("\t")[4] for line in lines[1:]]


def _run_15a(
    *,
    cohort_name: str,
    calls: list[list[int]],
    variants: list[tuple[str, int, str, str]],
    samples: list[str],
    covariates_rows: list[tuple[str, int, str]],
    tmp_path: Path,
) -> Path:
    cohort_dir = tmp_path / cohort_name
    cohort_dir.mkdir()
    bfile = cohort_dir / "geno"
    write_synthetic_plink(bfile, calls=calls, variants=variants, samples=samples)
    cov_path = cohort_dir / "covariates.txt"
    _write_covariates_file(cov_path, covariates_rows)

    out = cohort_dir / "out"
    log = cohort_dir / "log.txt"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--study-name", cohort_name,
            "--bfile", str(bfile),
            "--covariates", str(cov_path),
            "--output-dir", str(out),
            "--log-file", str(log),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    return out


def test_two_cohort_d_and_b_aggregation_match_stacked_reference(tmp_path):
    out_a = _run_15a(
        cohort_name="cohort_a",
        calls=COHORT_A_CALLS,
        variants=SHARED_VARIANTS,
        samples=COHORT_A_SAMPLES,
        covariates_rows=COHORT_A_COVARIATES,
        tmp_path=tmp_path,
    )
    out_b = _run_15a(
        cohort_name="cohort_b",
        calls=COHORT_B_CALLS,
        variants=SHARED_VARIANTS,
        samples=COHORT_B_SAMPLES,
        covariates_rows=COHORT_B_COVARIATES,
        tmp_path=tmp_path,
    )

    assert _read_variant_ids(out_a) == _read_variant_ids(out_b)

    d_a = np.load(out_a / "D.npy", allow_pickle=False)
    d_b = np.load(out_b / "D.npy", allow_pickle=False)
    b_a = np.load(out_a / "B.npy", allow_pickle=False)
    b_b = np.load(out_b / "B.npy", allow_pickle=False)

    x_a = np.asarray(COHORT_A_CALLS, dtype=np.float64).T
    x_b = np.asarray(COHORT_B_CALLS, dtype=np.float64).T
    c_a = _covariate_matrix(COHORT_A_COVARIATES)
    c_b = _covariate_matrix(COHORT_B_COVARIATES)

    x_full = np.vstack([x_a, x_b])
    c_full = np.vstack([c_a, c_b])

    d_ref = c_full.T @ c_full
    b_ref = x_full.T @ c_full

    np.testing.assert_allclose(d_a + d_b, d_ref, rtol=1e-12)
    np.testing.assert_allclose(b_a + b_b, b_ref, rtol=1e-12)


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

    # Intercept-only adjustment must reduce exactly to the Pearson correlation of
    # the stacked dosage matrix (the GoDMC mQTL estimand). Confirm the reference
    # equals np.corrcoef before masking to the LD window.
    np.testing.assert_allclose(R_ref, np.corrcoef(x_full, rowvar=False), rtol=1e-12)

    # The panel only stores pairs within the LD radius (DEFAULT_LD_RADIUS_BP =
    # 1,000,000 bp): pos_j <= pos_i + radius_bp. Variants here are at 1/2/3 Mbp,
    # so the (1Mbp, 3Mbp) pair is 2 Mbp apart and legitimately excluded. Mask the
    # dense reference to the same window before comparing the upper triangle.
    pos = np.array([v[1] for v in SHARED_VARIANTS])
    window = pos[None, :] <= pos[:, None] + 1_000_000
    R_ref = np.where(np.triu(window), R_ref, 0.0)

    R_got = _read_r_full(panel, "22", len(SHARED_VARIANTS))
    np.testing.assert_allclose(np.triu(R_got), np.triu(R_ref), rtol=1e-6, atol=1e-6)


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