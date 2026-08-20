"""End-to-end section-15 shell round-trip: ld_encrypt_cohort.sh (per-chromosome
packing) -> ld_decrypt_cohort.sh (decrypt + ld_reassemble_cohort.py merge).

Exercises the real shell entry points with real gpg, complementing
test_ld_reassemble_cohort.py (which drives the reassembler Python module directly
with Python-built archives). Fixtures mirror that module test's per-chromosome
layout so the merged output survives checksum verification.
"""

import gzip
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
ENCRYPT = REPO_ROOT / "resources" / "genetics" / "ld_encrypt_cohort.sh"
DECRYPT = REPO_ROOT / "resources" / "genetics" / "ld_decrypt_cohort.sh"
PASSPHRASE = "testpass"

sys.path.insert(0, str(REPO_ROOT / "resources" / "genetics"))
import ld_checksums as ck  # noqa: E402
from ld_qc import ordered_variant_digest  # noqa: E402


def _make_gpg_wrapper(tmp_path):
    """GPG wrapper with an isolated keyring."""
    gnupg_home = tmp_path / "gnupg"
    gnupg_home.mkdir(mode=0o700)
    wrapper = tmp_path / "gpg_batch.sh"
    wrapper.write_text(
        "#!/usr/bin/env bash\n"
        f'exec gpg --homedir "{gnupg_home}" "$@"\n'
    )
    wrapper.chmod(0o755)
    return wrapper


def _env(tmp_path, passphrase=PASSPHRASE):
    tmp_path.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["GPG"] = str(_make_gpg_wrapper(tmp_path))
    passphrase_file = tmp_path / "passphrase.txt"
    passphrase_file.write_text(passphrase + "\n")
    env["LD_GPG_PASSPHRASE_FILE"] = str(passphrase_file)
    return env


def _manifest(study, chrom, pos):
    variant_id = f"{chrom}:{pos}:A:C"
    digest = ordered_variant_digest([variant_id])
    return {
        "study_name": study,
        "generated_at_utc": "2026-01-01T00:00:00+00:00",
        "module": "15a",
        "status": "scaffold",
        "genome_build": "GRCh37",
        "autosomes_only": True,
        "input_bfile": "data",
        "sample_alignment": {"counts": {"final_sample_count": 10}},
        "variant_index": {
            "schema_version": "v0.4-canonical-row-alignment",
            "columns": ["chr", "pos", "ref", "alt", "variant_id",
                        "n_nonmissing", "n_imputed", "genotype_mean"],
            "chromosome_filter": chrom,
            "filters_applied": [f"chromosome {chrom}"],
            "counts": {
                "total_rows": 2,
                "kept_count": 1,
                "excluded_non_autosomal": 0,
                "excluded_other_chromosome": 1,
                "excluded_non_biallelic_snp": 0,
                "excluded_mhc_region": 0,
                "excluded_multiallelic_position": 0,
            },
            "n_samples_used": 10,
            "variant_id_digest": digest,
        },
        "covariate_schema": {
            "schema_id": "intercept-only-v1",
            "required_columns": [],
            "matrix_columns": ["intercept"],
            "sex_factor_recode": {},
        },
        "hail": {},
        "B_block": {"filename": "B.npy", "shape": [1, 1], "dtype": "float64",
                    "variant_id_digest": digest},
        "A_blocks": {
            "radius_bp": 1000000,
            "block_size": 4096,
            "chunk_rows": 50000,
            "max_dense_gb": 1.0,
            "variant_id_digest": digest,
            "chromosomes": {chrom: {
                "n_variants": 1,
                "n_chunks": 1,
                "variant_id_digest": digest,
                "chunks": [{"row_start": 0, "row_stop": 1,
                            "column_start": 0, "column_stop": 1,
                            "n_rows": 1, "n_cols": 1,
                            "row_variant_id_digest": digest,
                            "column_variant_id_digest": digest,
                            "directory": f"A_blocks/chr{chrom}/chunk_000000"}],
            }},
        },
        "notes": [],
    }


def _build_chromosome(cs, study, chrom, pos, b_value):
    """A single-chromosome 15a-style cohort_stats dir with real checksums."""
    a_chunk = cs / "A_blocks" / f"chr{chrom}" / "chunk_000000"
    a_chunk.mkdir(parents=True)
    (a_chunk / "values.txt").write_text(f"chr{chrom}\n", encoding="utf-8")
    np.save(cs / "B.npy", np.array([[b_value]], dtype=np.float64), allow_pickle=False)
    np.save(cs / "D.npy", np.array([[10.0]], dtype=np.float64), allow_pickle=False)
    with gzip.open(cs / "variants.tsv.gz", "wt") as handle:
        handle.write("chr\tpos\tref\talt\tvariant_id\tn_nonmissing\tn_imputed\tgenotype_mean\n")
        handle.write(f"{chrom}\t{pos}\tA\tC\t{chrom}:{pos}:A:C\t10\t0\t{b_value}\n")
    (cs / "manifest.json").write_text(
        json.dumps(_manifest(study, chrom, pos), indent=2) + "\n", encoding="utf-8")
    (cs / "qc_report.txt").write_text("ok\n", encoding="utf-8")
    ck.write_cohort_checksums(cs)
    return cs


def _encrypt(env, cohort, upload, study):
    return subprocess.run([str(ENCRYPT), str(cohort), str(upload), study],
                          env=env, capture_output=True, text=True)


def _decrypt(env, upload, merged, study, passphrase_file=None):
    passphrase_file = passphrase_file or env["LD_GPG_PASSPHRASE_FILE"]
    return subprocess.run(
        [str(DECRYPT), str(upload), str(merged), study, str(passphrase_file)],
        env=env,
        capture_output=True,
        text=True,
    )


def test_single_chromosome_roundtrip_produces_valid_cohort(tmp_path):
    env = _env(tmp_path)
    cohort = _build_chromosome(tmp_path / "chr22_cs", "cohortA", "22", 2201, 0.2)
    upload = tmp_path / "upload"
    merged = tmp_path / "merged"
    assert _encrypt(env, cohort, upload, "cohortA").returncode == 0
    r = _decrypt(env, upload, merged, "cohortA")
    assert r.returncode == 0, r.stderr

    # merged output is a valid, checksum-verified reassembled cohort
    ck.verify_cohort_checksums(merged)
    manifest = json.loads((merged / "manifest.json").read_text())
    assert manifest["study_name"] == "cohortA"
    assert manifest["status"] == "reassembled"
    assert (merged / "A_blocks" / "chr22" / "chunk_000000" / "values.txt").is_file()
    np.testing.assert_array_equal(np.load(merged / "B.npy"), np.array([[0.2]]))
    # input .aes archives are left untouched
    assert (upload / "cohortA_chr22_15_scaffold.tgz.aes").is_file()


def test_two_chromosome_roundtrip_merges_into_one_cohort(tmp_path):
    env = _env(tmp_path)
    upload = tmp_path / "upload"
    merged = tmp_path / "merged"
    # two per-chromosome cohort dirs, encrypted into the same upload area
    chr1 = _build_chromosome(tmp_path / "chr1_cs", "cohortA", "1", 101, 0.1)
    chr22 = _build_chromosome(tmp_path / "chr22_cs", "cohortA", "22", 2201, 0.2)
    assert _encrypt(env, chr1, upload, "cohortA").returncode == 0
    assert _encrypt(env, chr22, upload, "cohortA").returncode == 0

    r = _decrypt(env, upload, merged, "cohortA")
    assert r.returncode == 0, r.stderr

    ck.verify_cohort_checksums(merged)
    manifest = json.loads((merged / "manifest.json").read_text())
    assert list(manifest["A_blocks"]["chromosomes"]) == ["1", "22"]
    assert manifest["variant_index"]["counts"]["kept_count"] == 2
    assert manifest["B_block"]["shape"] == [2, 1]
    assert (merged / "A_blocks" / "chr1" / "chunk_000000" / "values.txt").is_file()
    assert (merged / "A_blocks" / "chr22" / "chunk_000000" / "values.txt").is_file()
    np.testing.assert_array_equal(np.load(merged / "B.npy"), np.array([[0.1], [0.2]]))


def test_each_cohort_uses_its_own_passphrase_file(tmp_path):
    env_a = _env(tmp_path / "gpg_a", "cohort-a-passphrase")
    env_b = _env(tmp_path / "gpg_b", "cohort-b-passphrase")
    upload = tmp_path / "upload"
    cohort_a = _build_chromosome(
        tmp_path / "cohort_a", "cohortA", "22", 2201, 0.2
    )
    cohort_b = _build_chromosome(
        tmp_path / "cohort_b", "cohortB", "22", 2201, 0.3
    )

    assert _encrypt(env_a, cohort_a, upload, "cohortA").returncode == 0
    assert _encrypt(env_b, cohort_b, upload, "cohortB").returncode == 0

    merged_a = tmp_path / "merged_a"
    merged_b = tmp_path / "merged_b"
    assert _decrypt(env_a, upload, merged_a, "cohortA").returncode == 0
    assert _decrypt(env_b, upload, merged_b, "cohortB").returncode == 0
    manifest_a = json.loads((merged_a / "manifest.json").read_text())
    manifest_b = json.loads((merged_b / "manifest.json").read_text())
    assert manifest_a["study_name"] == "cohortA"
    assert manifest_b["study_name"] == "cohortB"

    wrong = _decrypt(
        env_a,
        upload,
        tmp_path / "wrong_passphrase",
        "cohortA",
        env_b["LD_GPG_PASSPHRASE_FILE"],
    )
    assert wrong.returncode != 0
    assert not (tmp_path / "wrong_passphrase" / "manifest.json").exists()


def test_md5_mismatch_fails_and_does_not_publish(tmp_path):
    env = _env(tmp_path)
    cohort = _build_chromosome(tmp_path / "chr22_cs", "cohortA", "22", 2201, 0.2)
    upload = tmp_path / "upload"
    merged = tmp_path / "merged"
    assert _encrypt(env, cohort, upload, "cohortA").returncode == 0
    # Corrupt the recorded plaintext md5 for the chunk so verification fails.
    md5 = upload / "cohortA_chr22_15_chr22_chunk_000000.md5sum"
    md5.write_text("0" * 32 + "  cohortA_chr22_15_chr22_chunk_000000.tgz\n")
    r = _decrypt(env, upload, merged, "cohortA")
    assert r.returncode != 0
    # nothing published to the merged output on failure
    assert not (merged / "manifest.json").exists()


def test_wrong_study_name_finds_no_archives(tmp_path):
    env = _env(tmp_path)
    cohort = _build_chromosome(tmp_path / "chr22_cs", "cohortA", "22", 2201, 0.2)
    upload = tmp_path / "upload"
    merged = tmp_path / "merged"
    assert _encrypt(env, cohort, upload, "cohortA").returncode == 0
    r = _decrypt(env, upload, merged, "cohortB")
    assert r.returncode != 0
    assert "scaffold" in r.stderr.lower()
