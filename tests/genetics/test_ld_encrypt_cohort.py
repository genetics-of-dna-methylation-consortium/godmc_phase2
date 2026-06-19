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
    return wrapper, gnupg_home


def _make_cohort(tmp_path, chunks=(("chr1", "chunk_0"), ("chr1", "chunk_1"), ("chr2", "chunk_0"))):
    cs = tmp_path / "cohort_stats"
    (cs / "A_blocks").mkdir(parents=True)
    (cs / "manifest.json").write_text('{"module": "15a"}')
    with gzip.open(cs / "variants.tsv.gz", "wt") as fh:
        fh.write("chr\tpos\n1\t1000\n")
    (cs / "D.npy").write_bytes(b"D-matrix-bytes")
    (cs / "B.npy").write_bytes(b"B-matrix-bytes")
    for chrom, chunk in chunks:
        d = cs / "A_blocks" / chrom / chunk
        d.mkdir(parents=True)
        (d / "part-00000").write_bytes(f"{chrom}/{chunk}/data".encode())
        (d / "metadata.json").write_text('{"block":1}')
    return cs


def _run(helper_env, cohort, out, study="testcohort"):
    return subprocess.run(
        [str(HELPER), str(cohort), str(out), study],
        env=helper_env, capture_output=True, text=True,
    )


def _gpg_env(tmp_path):
    wrapper, _ = _make_gpg_wrapper(tmp_path)
    env = dict(os.environ)
    env["GPG"] = str(wrapper)
    return env, wrapper


def test_encrypts_scaffold_and_each_chunk(tmp_path):
    cohort = _make_cohort(tmp_path)
    out = tmp_path / "upload"
    env, _ = _gpg_env(tmp_path)
    r = _run(env, cohort, out)
    assert r.returncode == 0, r.stderr
    # scaffold
    assert (out / "testcohort_15_scaffold.tgz.aes").is_file()
    assert (out / "testcohort_15_scaffold.md5sum").is_file()
    # one archive per chunk
    for base in ("testcohort_15_chr1_chunk_0", "testcohort_15_chr1_chunk_1",
                 "testcohort_15_chr2_chunk_0"):
        assert (out / f"{base}.tgz.aes").is_file(), base
        assert (out / f"{base}.md5sum").is_file(), base
    # no plaintext tarballs left behind
    assert list(out.glob("*.tgz")) == []


def test_chunk_roundtrip_reproduces_tree(tmp_path):
    cohort = _make_cohort(tmp_path)
    out = tmp_path / "upload"
    wrapper, _ = _make_gpg_wrapper(tmp_path)
    env = dict(os.environ); env["GPG"] = str(wrapper)
    assert _run(env, cohort, out).returncode == 0
    aes = out / "testcohort_15_chr1_chunk_0.tgz.aes"
    dec = out / "decrypted.tgz"
    subprocess.run([str(wrapper), "--output", str(dec), "-d", str(aes)], check=True)
    extract = tmp_path / "extract"; extract.mkdir()
    with tarfile.open(dec) as tf:
        tf.extractall(extract)
    assert (extract / "chr1" / "chunk_0" / "part-00000").read_bytes() == b"chr1/chunk_0/data"
