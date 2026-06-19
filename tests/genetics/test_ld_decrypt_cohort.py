import gzip
import hashlib
import os
import subprocess
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
ENCRYPT = REPO_ROOT / "resources" / "genetics" / "ld_encrypt_cohort.sh"
DECRYPT = REPO_ROOT / "resources" / "genetics" / "ld_decrypt_cohort.sh"
PASSPHRASE = "testpass"


def _make_gpg_wrapper(tmp_path):
    """Non-interactive gpg wrapper (loopback passphrase) + isolated keyring."""
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


def _make_cohort(tmp_path, study="testcohort",
                 chunks=(("chr1", "chunk_0"), ("chr1", "chunk_1"), ("chr2", "chunk_0"))):
    cs = tmp_path / "cohort_stats"
    (cs / "A_blocks").mkdir(parents=True)
    (cs / "manifest.json").write_text(f'{{"study_name": "{study}"}}')
    with gzip.open(cs / "variants.tsv.gz", "wt") as fh:
        fh.write("chr\tpos\n1\t1000\n")
    (cs / "D.npy").write_bytes(b"D-matrix-bytes")
    (cs / "B.npy").write_bytes(b"B-matrix-bytes")
    (cs / "checksums.json").write_text('{"algorithm": "blake2b", "files": {}}')
    for chrom, chunk in chunks:
        d = cs / "A_blocks" / chrom / chunk
        d.mkdir(parents=True)
        (d / "part-00000").write_bytes(f"{chrom}/{chunk}/data".encode())
        (d / "metadata.json").write_text('{"block":1}')
    return cs


def _digests(root):
    """Map every file's POSIX relpath under root -> sha256 hex (for tree compare)."""
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            out[p.relative_to(root).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def _env(wrapper):
    env = dict(os.environ)
    env["GPG"] = str(wrapper)
    return env


def _encrypt(env, cohort, upload, study="testcohort"):
    return subprocess.run([str(ENCRYPT), str(cohort), str(upload), study],
                          env=env, capture_output=True, text=True)


def _decrypt(env, upload, rebuilt, study="testcohort"):
    return subprocess.run([str(DECRYPT), str(upload), str(rebuilt), study],
                          env=env, capture_output=True, text=True)


def test_roundtrip_reproduces_cohort_tree(tmp_path):
    wrapper = _make_gpg_wrapper(tmp_path)
    env = _env(wrapper)
    cohort = _make_cohort(tmp_path)
    upload = tmp_path / "upload"
    rebuilt = tmp_path / "rebuilt"
    assert _encrypt(env, cohort, upload).returncode == 0
    r = _decrypt(env, upload, rebuilt)
    assert r.returncode == 0, r.stderr
    # rebuilt tree matches the original byte-for-byte (ignoring the .staging workdir)
    rebuilt_digests = {k: v for k, v in _digests(rebuilt).items()
                       if not k.startswith(".staging/")}
    assert rebuilt_digests == _digests(cohort)
    # no leftover plaintext tarballs anywhere under rebuilt
    assert list(rebuilt.rglob("*.tgz")) == []
    # input .aes archives left untouched
    assert (upload / "testcohort_15_scaffold.tgz.aes").is_file()


def test_resume_skips_already_restored_chunk(tmp_path):
    wrapper = _make_gpg_wrapper(tmp_path)
    env = _env(wrapper)
    cohort = _make_cohort(tmp_path)
    upload = tmp_path / "upload"
    rebuilt = tmp_path / "rebuilt"
    assert _encrypt(env, cohort, upload).returncode == 0
    assert _decrypt(env, upload, rebuilt).returncode == 0
    chunk = rebuilt / "A_blocks" / "chr1" / "chunk_0" / "part-00000"
    mtime_before = chunk.stat().st_mtime_ns
    r2 = _decrypt(env, upload, rebuilt)
    assert r2.returncode == 0, r2.stderr
    assert "skip testcohort_15_chr1_chunk_0" in r2.stdout
    assert chunk.stat().st_mtime_ns == mtime_before  # not re-extracted


def test_md5_mismatch_fails_before_untar(tmp_path):
    wrapper = _make_gpg_wrapper(tmp_path)
    env = _env(wrapper)
    cohort = _make_cohort(tmp_path)
    upload = tmp_path / "upload"
    rebuilt = tmp_path / "rebuilt"
    assert _encrypt(env, cohort, upload).returncode == 0
    # Corrupt the recorded md5 for one chunk so md5sum -c fails.
    md5 = upload / "testcohort_15_chr1_chunk_0.md5sum"
    md5.write_text("0" * 32 + "  testcohort_15_chr1_chunk_0.tgz\n")
    r = _decrypt(env, upload, rebuilt)
    assert r.returncode != 0
    # the corrupted chunk must NOT have been published into the tree
    assert not (rebuilt / "A_blocks" / "chr1" / "chunk_0").exists()


def test_study_name_mismatch_fails(tmp_path):
    wrapper = _make_gpg_wrapper(tmp_path)
    env = _env(wrapper)
    cohort = _make_cohort(tmp_path, study="cohortA")
    upload = tmp_path / "upload"
    rebuilt = tmp_path / "rebuilt"
    assert _encrypt(env, cohort, upload, study="cohortA").returncode == 0
    # Decrypt asking for a different study_name than the manifest records.
    r = subprocess.run([str(DECRYPT), str(upload), str(rebuilt), "cohortB"],
                       env=env, capture_output=True, text=True)
    assert r.returncode == 1
    assert "study_name" in r.stderr and "cohortA" in r.stderr


def test_zero_chunks_fails(tmp_path):
    wrapper = _make_gpg_wrapper(tmp_path)
    env = _env(wrapper)
    cohort = _make_cohort(tmp_path, chunks=())  # scaffold only, no chunks
    upload = tmp_path / "upload"
    rebuilt = tmp_path / "rebuilt"
    # encrypt fails on zero chunks, so stage only the scaffold archive by hand:
    upload.mkdir()
    subprocess.run([str(ENCRYPT), str(cohort), str(upload), "testcohort"],
                   env=env, capture_output=True, text=True)  # produces scaffold, then errors
    r = _decrypt(env, upload, rebuilt)
    assert r.returncode == 1
    assert "no A_blocks chunk archives" in r.stderr
