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
