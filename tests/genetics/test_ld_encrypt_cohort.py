import gzip
import os
import subprocess
import tarfile
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
HELPER = REPO_ROOT / "resources" / "genetics" / "ld_encrypt_cohort.sh"
PASSPHRASE = "testpass"


def _make_gpg_wrapper(tmp_path):
    """A non-interactive gpg wrapper (loopback passphrase) + isolated keyring."""
    gnupg_home = tmp_path / "gnupg"
    gnupg_home.mkdir(mode=0o700)
    wrapper = tmp_path / "gpg_batch.sh"
    wrapper.write_text(
        "#!/usr/bin/env bash\n"
        f'exec gpg --homedir "{gnupg_home}" --batch --yes '
        f'--pinentry-mode loopback --passphrase "{PASSPHRASE}" "$@"\n'
    )
    wrapper.chmod(0o755)
    return wrapper


def _make_chromosome(tmp_path, chrom="22", chunks=("chunk_0", "chunk_1")):
    """A single-chromosome 15a-style cohort dir (one A_blocks/chr<C> subtree)."""
    cs = tmp_path / "cohort_stats"
    a = cs / "A_blocks" / f"chr{chrom}"
    a.mkdir(parents=True)
    (cs / "manifest.json").write_text('{"module": "15a"}')
    with gzip.open(cs / "variants.tsv.gz", "wt") as fh:
        fh.write("chr\tpos\n{}\t1000\n".format(chrom))
    (cs / "D.npy").write_bytes(b"D-matrix-bytes")
    (cs / "B.npy").write_bytes(b"B-matrix-bytes")
    (cs / "checksums.json").write_text('{"algorithm": "blake2b", "files": {}}')
    (cs / "qc_report.txt").write_text("ok\n")
    for chunk in chunks:
        d = a / chunk
        d.mkdir(parents=True)
        (d / "part-00000").write_bytes(f"chr{chrom}/{chunk}/data".encode())
        (d / "metadata.json").write_text('{"block":1}')
    return cs


def _run(env, cohort, out, study="testcohort"):
    return subprocess.run(
        [str(HELPER), str(cohort), str(out), study],
        env=env, capture_output=True, text=True,
    )


def _env(tmp_path):
    env = dict(os.environ)
    env["GPG"] = str(_make_gpg_wrapper(tmp_path))
    passphrase_file = tmp_path / "passphrase.txt"
    passphrase_file.write_text(PASSPHRASE + "\n")
    env["LD_GPG_PASSPHRASE_FILE"] = str(passphrase_file)
    return env, Path(env["GPG"])


def test_parses():
    r = subprocess.run(["bash", "-n", str(HELPER)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_emits_per_chromosome_scaffold_and_chunks(tmp_path):
    cohort = _make_chromosome(tmp_path, chrom="22")
    out = tmp_path / "upload"
    env, _ = _env(tmp_path)
    r = _run(env, cohort, out)
    assert r.returncode == 0, r.stderr
    # scaffold carries the chromosome in its name (reassembler contract)
    assert (out / "testcohort_chr22_15_scaffold.tgz.aes").is_file()
    assert (out / "testcohort_chr22_15_scaffold.md5sum").is_file()
    # one archive per chunk, per-chromosome naming
    for base in ("testcohort_chr22_15_chr22_chunk_0",
                 "testcohort_chr22_15_chr22_chunk_1"):
        assert (out / f"{base}.tgz.aes").is_file(), base
        assert (out / f"{base}.md5sum").is_file(), base
    # no plaintext tarballs left behind
    assert list(out.glob("*.tgz")) == []


def test_scaffold_archive_includes_checksums_and_qc(tmp_path):
    cohort = _make_chromosome(tmp_path, chrom="22")
    out = tmp_path / "upload"
    env, gpg = _env(tmp_path)
    assert _run(env, cohort, out).returncode == 0
    dec = out / "scaffold.tgz"
    subprocess.run([str(gpg), "--output", str(dec), "-d",
                    str(out / "testcohort_chr22_15_scaffold.tgz.aes")], check=True)
    with tarfile.open(dec) as tf:
        names = set(tf.getnames())
    assert {"manifest.json", "variants.tsv.gz", "D.npy", "B.npy",
            "checksums.json", "qc_report.txt"} <= names


def test_chunk_roundtrip_reproduces_tree(tmp_path):
    cohort = _make_chromosome(tmp_path, chrom="22")
    out = tmp_path / "upload"
    env, gpg = _env(tmp_path)
    assert _run(env, cohort, out).returncode == 0
    dec = out / "decrypted.tgz"
    subprocess.run([str(gpg), "--output", str(dec), "-d",
                    str(out / "testcohort_chr22_15_chr22_chunk_0.tgz.aes")], check=True)
    extract = tmp_path / "extract"
    extract.mkdir()
    with tarfile.open(dec) as tf:
        tf.extractall(extract)
    assert (extract / "chr22" / "chunk_0" / "part-00000").read_bytes() == b"chr22/chunk_0/data"


def test_resume_skips_already_encrypted(tmp_path):
    cohort = _make_chromosome(tmp_path, chrom="22")
    out = tmp_path / "upload"
    env, _ = _env(tmp_path)
    assert _run(env, cohort, out).returncode == 0
    aes = out / "testcohort_chr22_15_chr22_chunk_0.tgz.aes"
    mtime_before = aes.stat().st_mtime_ns
    r2 = _run(env, cohort, out)
    assert r2.returncode == 0, r2.stderr
    assert "skip" in r2.stdout.lower()
    assert aes.stat().st_mtime_ns == mtime_before  # not regenerated


def test_zero_chunks_fails(tmp_path):
    cohort = _make_chromosome(tmp_path, chrom="22", chunks=())
    out = tmp_path / "upload"
    env, _ = _env(tmp_path)
    r = _run(env, cohort, out)
    assert r.returncode != 0
    assert "no A_blocks chunks" in r.stderr


def test_multiple_chromosomes_fails(tmp_path):
    cohort = _make_chromosome(tmp_path, chrom="22")
    # add a second chromosome subtree -> ambiguous single-chromosome input
    extra = cohort / "A_blocks" / "chr1" / "chunk_0"
    extra.mkdir(parents=True)
    (extra / "part-00000").write_bytes(b"chr1/chunk_0/data")
    out = tmp_path / "upload"
    env, _ = _env(tmp_path)
    r = _run(env, cohort, out)
    assert r.returncode != 0
    assert "one chromosome" in r.stderr.lower()


def test_manifest_chunk_count_mismatch_warns_but_succeeds(tmp_path):
    if not _have_jq():
        pytest.skip("jq not available")
    cohort = _make_chromosome(tmp_path, chrom="22")  # 2 chunks on disk
    (cohort / "manifest.json").write_text(
        '{"A_blocks": {"chromosomes": {"22": {"n_chunks": 5}}}}'
    )
    out = tmp_path / "upload"
    env, _ = _env(tmp_path)
    r = _run(env, cohort, out)
    assert r.returncode == 0, r.stderr
    assert "WARNING" in r.stderr and "5" in r.stderr and "2" in r.stderr


def _have_jq():
    return subprocess.run(["bash", "-c", "command -v jq"],
                          capture_output=True).returncode == 0
