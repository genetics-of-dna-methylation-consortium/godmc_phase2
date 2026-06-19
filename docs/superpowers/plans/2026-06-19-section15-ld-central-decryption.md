# Section 15 (LD) Central Decryption Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a central-side `ld_decrypt_cohort.sh` that decrypts a cohort's uploaded section-15 GPG archives, verifies integrity, and rebuilds the `cohort_stats/` tree that `ld_aggregate.py accumulate` consumes.

**Architecture:** A standalone bash script that mirrors `ld_encrypt_cohort.sh` (no network I/O). It decrypts each `.aes`, verifies the sibling `.md5sum` against the plaintext tar, then atomically untars into the rebuilt tree. One cohort-side edit adds `checksums.json` to the scaffold tar so the blake2b integrity layer survives the round-trip and `accumulate` can re-verify it.

**Tech Stack:** Bash (`set -euo pipefail`), GnuPG symmetric AES256, `tar`/`gzip`, `md5sum`; pytest with a real-gpg loopback wrapper for tests; Python (`ld_aggregate`, `ld_checksums`) for the end-to-end test.

## Global Constraints

- Decryption uses symmetric `gpg`; the binary is overridable via the `GPG` env var (default `gpg`) — copied verbatim from `ld_encrypt_cohort.sh`.
- No network I/O in the decrypt script (mirror of the encrypt script).
- The rebuilt `output_dir` must equal what `accumulate --cohort-dir` consumes: top-level `manifest.json`, `variants.tsv.gz`, `B.npy`, `D.npy`, `checksums.json`, and `A_blocks/chr<C>/chunk_<N>/`.
- `study_name` is NOT an `accumulate` argument — `accumulate` reads it from `manifest.json["study_name"]` (`resources/genetics/ld_aggregate.py:440`). In decrypt it is only an archive-selection prefix plus a manifest cross-check.
- Both integrity layers are preserved: `.md5sum` (md5 of plaintext tar, checked at decrypt) and `checksums.json` (blake2b, checked by `accumulate`).
- Spec: `docs/superpowers/specs/2026-06-19-section15-ld-central-decryption-design.md`.

## File Structure

- `resources/genetics/ld_encrypt_cohort.sh` — **modify** (Task 1): add `checksums.json` to the scaffold tar (line 47).
- `resources/genetics/ld_decrypt_cohort.sh` — **create** (Task 2): the decrypt + reassemble script.
- `tests/genetics/test_ld_encrypt_cohort.py` — **modify** (Task 1): fixture writes `checksums.json`; assert it is in the scaffold.
- `tests/genetics/test_ld_decrypt_cohort.py` — **create** (Tasks 2 & 3): round-trip, resume, md5-fail, study-name-mismatch, zero-chunks, and end-to-end-with-accumulate tests.

---

## Task 1: Add `checksums.json` to the scaffold archive (cohort-side)

**Files:**
- Modify: `resources/genetics/ld_encrypt_cohort.sh:46-47`
- Test: `tests/genetics/test_ld_encrypt_cohort.py` (fixture + new assertion)

**Interfaces:**
- Consumes: nothing new. `checksums.json` is written by 15a's final step (`ld_prepare_stats.py:328` → `ld_checksums.write_cohort_checksums`), so it is present in `cohort_stats_dir` at encrypt time.
- Produces: the scaffold archive `<study>_15_scaffold.tgz.aes` now contains `checksums.json` alongside `manifest.json variants.tsv.gz D.npy B.npy`. Task 2 and Task 3 rely on this.

- [ ] **Step 1: Update the encrypt-test fixture to include `checksums.json`**

In `tests/genetics/test_ld_encrypt_cohort.py`, edit `_make_cohort` to write a `checksums.json` file (a placeholder is fine here — Task 1 only checks packaging, not blake2b validity). Add this line just before `return cs` (after the `B.npy` write at line 35):

```python
    (cs / "checksums.json").write_text('{"algorithm": "blake2b", "files": {}}')
```

- [ ] **Step 2: Write the failing assertion that the scaffold contains `checksums.json`**

Add this test to `tests/genetics/test_ld_encrypt_cohort.py`:

```python
def test_scaffold_archive_includes_checksums_json(tmp_path):
    cohort = _make_cohort(tmp_path)
    out = tmp_path / "upload"
    wrapper, _ = _make_gpg_wrapper(tmp_path)
    env = dict(os.environ); env["GPG"] = str(wrapper)
    assert _run(env, cohort, out).returncode == 0
    dec = out / "scaffold.tgz"
    subprocess.run(
        [str(wrapper), "--output", str(dec), "-d",
         str(out / "testcohort_15_scaffold.tgz.aes")], check=True)
    with tarfile.open(dec) as tf:
        names = tf.getnames()
    assert "checksums.json" in names, names
    assert "manifest.json" in names
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `pytest tests/genetics/test_ld_encrypt_cohort.py::test_scaffold_archive_includes_checksums_json -v`
Expected: FAIL — `assert "checksums.json" in names` (scaffold tar currently packs only `manifest.json variants.tsv.gz D.npy B.npy`).

- [ ] **Step 4: Add `checksums.json` to the scaffold tar**

In `resources/genetics/ld_encrypt_cohort.sh`, change the scaffold `tar` (lines 46-47) from:

```bash
  tar czf "${output_dir}/${scaffold}.tgz" -C "${cohort_stats_dir}" \
    manifest.json variants.tsv.gz D.npy B.npy
```

to:

```bash
  tar czf "${output_dir}/${scaffold}.tgz" -C "${cohort_stats_dir}" \
    manifest.json variants.tsv.gz D.npy B.npy checksums.json
```

- [ ] **Step 5: Run the new test and the full encrypt suite to verify they pass**

Run: `pytest tests/genetics/test_ld_encrypt_cohort.py -v`
Expected: PASS — all existing tests plus `test_scaffold_archive_includes_checksums_json`.

- [ ] **Step 6: Commit**

```bash
git add resources/genetics/ld_encrypt_cohort.sh tests/genetics/test_ld_encrypt_cohort.py
git commit -m "feat(15): include checksums.json in scaffold archive for round-trip integrity"
```

---

## Task 2: Create `ld_decrypt_cohort.sh` (decrypt + reassemble)

**Files:**
- Create: `resources/genetics/ld_decrypt_cohort.sh`
- Test: `tests/genetics/test_ld_decrypt_cohort.py`

**Interfaces:**
- Consumes: the `.aes`/`.md5sum` archive set produced by `ld_encrypt_cohort.sh` (scaffold now including `checksums.json` from Task 1). Archive naming: `<study>_15_scaffold` and `<study>_15_chr<C>_chunk_<N>`; chunk tars are built with `-C <cohort_stats>/A_blocks` so their contents are `chr<C>/chunk_<N>/…`. Each `.md5sum` references the bare filename `<base>.tgz`.
- Produces: CLI `ld_decrypt_cohort.sh <input_dir> <output_dir> <study_name>` that rebuilds the `cohort_stats/` tree at `output_dir`. Exit codes: `0` success; `2` bad args; `1` zero chunks found or study-name mismatch. Honours the `GPG` env override.

- [ ] **Step 1: Write the failing round-trip test**

Create `tests/genetics/test_ld_decrypt_cohort.py`. It reuses the real-gpg loopback wrapper pattern from the encrypt test, encrypts a fixture cohort, then decrypts it and asserts a byte-exact tree:

```python
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
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pytest tests/genetics/test_ld_decrypt_cohort.py::test_roundtrip_reproduces_cohort_tree -v`
Expected: FAIL — `ld_decrypt_cohort.sh` does not exist yet (non-zero return / file-not-found).

- [ ] **Step 3: Create `ld_decrypt_cohort.sh` with the full restore logic**

Create `resources/genetics/ld_decrypt_cohort.sh` with exactly this content:

```bash
#!/usr/bin/env bash
# ld_decrypt_cohort.sh — decrypt + reassemble section-15 (LD) cohort archives
# into the cohort_stats/ tree that `ld_aggregate accumulate` consumes. No network
# I/O. Mirror of ld_encrypt_cohort.sh.
# See docs/superpowers/specs/2026-06-19-section15-ld-central-decryption-design.md
#
# Usage: ld_decrypt_cohort.sh <input_dir> <output_dir> <study_name>
#   input_dir   dir holding <study>_15_*.tgz.aes + .md5sum (left untouched)
#   output_dir  where the cohort_stats/ tree is rebuilt (== accumulate --cohort-dir)
#   study_name  selects this cohort's archives by filename prefix
#
# Override the gpg binary/wrapper for testing via the GPG env var.
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input_dir> <output_dir> <study_name>" >&2
  exit 2
fi

# Absolutise paths: md5sum -c runs from the staging dir, so a relative input_dir
# would otherwise break. output_dir is created first so it can be resolved too.
input_dir="$(cd "$1" && pwd)"
mkdir -p "$2"
output_dir="$(cd "$2" && pwd)"
study_name="$3"
GPG="${GPG:-gpg}"

staging="${output_dir}/.staging"
mkdir -p "${staging}"

# Decrypt ${input_dir}/$1.tgz.aes -> ${staging}/$1.tgz and verify its .md5sum.
# The input .aes is never modified or removed. The caller derives the tar path as
# "${staging}/$1.tgz" itself (do NOT capture this function's stdout — md5sum -c
# prints an "OK" line there).
decrypt_and_verify () {
  local base="$1"
  local tgz="${staging}/${base}.tgz"
  "${GPG}" --output "${tgz}" --decrypt "${input_dir}/${base}.tgz.aes"
  ( cd "${staging}" && md5sum -c "${input_dir}/${base}.md5sum" )
}

# --- 1. Scaffold (small files: manifest, variants, D, B, checksums.json) ---
scaffold="${study_name}_15_scaffold"
if [ -f "${output_dir}/manifest.json" ]; then
  echo "[ld_decrypt] skip ${scaffold} (manifest.json already present)"
else
  decrypt_and_verify "${scaffold}"
  tgz="${staging}/${scaffold}.tgz"
  tmp="${staging}/scaffold_extract.$$"
  rm -rf "${tmp}"; mkdir -p "${tmp}"
  tar xzf "${tgz}" -C "${tmp}"
  # Move every scaffold file into place, manifest.json LAST so its presence is a
  # reliable "scaffold complete" resume marker.
  for f in "${tmp}"/*; do
    [ "$(basename "${f}")" = "manifest.json" ] && continue
    mv -f "${f}" "${output_dir}/"
  done
  mv -f "${tmp}/manifest.json" "${output_dir}/manifest.json"
  rm -rf "${tmp}" "${tgz}"
  echo "[ld_decrypt] restored ${scaffold}"
fi

# --- 2. Study-name cross-check (hard fail on mismatch) ---
got_study="$(python -c 'import json,sys; print(json.load(open(sys.argv[1])).get("study_name",""))' \
  "${output_dir}/manifest.json")"
if [ "${got_study}" != "${study_name}" ]; then
  echo "[ld_decrypt] ERROR: manifest study_name '${got_study}' != requested '${study_name}'" >&2
  exit 1
fi

# --- 3. Per-chunk A_blocks archives ---
mkdir -p "${output_dir}/A_blocks"
n_chunks=0
shopt -s nullglob
for aes in "${input_dir}/${study_name}_15_chr"*"_chunk_"*".tgz.aes"; do
  base="$(basename "${aes}" .tgz.aes)"           # e.g. study_15_chr1_chunk_0
  suffix="${base#${study_name}_15_}"             # chr1_chunk_0
  chr_name="${suffix%%_chunk_*}"                 # chr1
  chunk_name="chunk_${suffix#*_chunk_}"          # chunk_0
  target="${output_dir}/A_blocks/${chr_name}/${chunk_name}"
  n_chunks=$((n_chunks + 1))
  if [ -d "${target}" ]; then
    echo "[ld_decrypt] skip ${base} (already restored)"
    continue
  fi
  decrypt_and_verify "${base}"
  tgz="${staging}/${base}.tgz"
  tmp="${staging}/chunk_extract.$$"
  rm -rf "${tmp}"; mkdir -p "${tmp}"
  tar xzf "${tgz}" -C "${tmp}"                    # tmp now holds chr<C>/chunk_<N>/...
  mkdir -p "${output_dir}/A_blocks/${chr_name}"
  mv "${tmp}/${chr_name}/${chunk_name}" "${target}"   # atomic publish
  rm -rf "${tmp}" "${tgz}"
  echo "[ld_decrypt] restored ${base}"
done
shopt -u nullglob

if [ "${n_chunks}" -eq 0 ]; then
  echo "[ld_decrypt] ERROR: no A_blocks chunk archives found for '${study_name}' in ${input_dir}" >&2
  exit 1
fi

# Best-effort integrity guard: compare against the manifest's declared chunk count.
manifest="${output_dir}/manifest.json"
if command -v jq >/dev/null 2>&1 && [ -f "${manifest}" ]; then
  expected="$(jq '[.A_blocks.chromosomes[].n_chunks] | add // 0' "${manifest}" 2>/dev/null || true)"
  if [ -n "${expected}" ] && [ "${expected}" -gt 0 ] && [ "${expected}" != "${n_chunks}" ]; then
    echo "[ld_decrypt] WARNING: manifest expects ${expected} chunks but restored ${n_chunks}" >&2
  fi
fi

rmdir "${staging}" 2>/dev/null || true
echo "[ld_decrypt] done: ${n_chunks} chunk archive(s) + scaffold in ${output_dir}"
```

Then make it executable:

```bash
chmod +x resources/genetics/ld_decrypt_cohort.sh
```

- [ ] **Step 4: Run the round-trip test to verify it passes**

Run: `pytest tests/genetics/test_ld_decrypt_cohort.py::test_roundtrip_reproduces_cohort_tree -v`
Expected: PASS — rebuilt tree byte-matches the original, no leftover `.tgz`, `.aes` untouched.

- [ ] **Step 5: Add the resume, md5-fail, study-name-mismatch, and zero-chunk tests**

Append to `tests/genetics/test_ld_decrypt_cohort.py`:

```python
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
```

- [ ] **Step 6: Run the full decrypt suite to verify it passes**

Run: `pytest tests/genetics/test_ld_decrypt_cohort.py -v`
Expected: PASS for all five tests.

- [ ] **Step 7: Commit**

```bash
git add resources/genetics/ld_decrypt_cohort.sh tests/genetics/test_ld_decrypt_cohort.py
git commit -m "feat(15): central LD decryption + reassembly (ld_decrypt_cohort.sh)"
```

---

## Task 3: End-to-end — decrypt feeds `accumulate` with blake2b verification on

**Files:**
- Test: `tests/genetics/test_ld_decrypt_cohort.py` (one new test)

**Interfaces:**
- Consumes: `ld_decrypt_cohort.sh` (Task 2), `ld_aggregate.accumulate(cohort_dir, precursor_dir, chunk_reader=..., verify_checksums=True)`, and `ld_checksums.write_cohort_checksums` (`resources/genetics/ld_checksums.py`).
- Produces: proof that a real, checksummed cohort survives encrypt → decrypt → `accumulate` with `verify_checksums=True` (the default) — i.e. the `checksums.json` round-trips and the blake2b layer validates against the rebuilt tree.

- [ ] **Step 1: Write the failing end-to-end test**

This builds a real (tiny) cohort with a valid `checksums.json` via `ld_checksums.write_cohort_checksums`, encrypts and decrypts it, then accumulates the rebuilt tree. It imports `ld_aggregate`/`ld_checksums` the same way `tests/genetics/test_ld_aggregate.py` does (the genetics test dir already resolves these imports). Append to `tests/genetics/test_ld_decrypt_cohort.py`:

```python
import json
import sys

import numpy as np

sys.path.insert(0, str(REPO_ROOT / "resources" / "genetics"))
import ld_aggregate as agg          # noqa: E402
import ld_checksums as ck           # noqa: E402


def _build_real_cohort(tmp_path, study="cohortA"):
    """A tiny 2-variant chr1 cohort with a real checksums.json + an A_blocks chunk."""
    cs = tmp_path / "cohort_stats"
    chunk = cs / "A_blocks" / "chr1" / "chunk_000000"
    chunk.mkdir(parents=True)
    (chunk / "part-00000").write_bytes(b"blockmatrix-bytes")
    rows = [("1", 100, "G", "A", "1:100:G:A", 4, 0, 1.0),
            ("1", 200, "C", "T", "1:200:C:T", 4, 0, 0.5)]
    with gzip.open(cs / "variants.tsv.gz", "wt") as fh:
        fh.write("chr\tpos\tref\talt\tvariant_id\tn_nonmissing\tn_imputed\tgenotype_mean\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    np.save(cs / "B.npy", np.array([[4.0], [2.0]]))
    np.save(cs / "D.npy", np.array([[4.0]]))
    manifest = {
        "study_name": study, "genome_build": "GRCh37",
        "covariate_schema": {"schema_id": "intercept_only",
            "matrix_columns": ["intercept"], "sex_factor_recode": {}},
        "variant_index": {"schema_version": "v0.3-with-genotype-stats"},
        "A_blocks": {"radius_bp": 1_000_000, "block_size": 8,
            "chromosomes": {"1": {"chunks": [
                {"row_start": 0, "row_stop": 2, "column_start": 0,
                 "column_stop": 2, "directory": "A_blocks/chr1/chunk_000000"}]}}},
    }
    (cs / "manifest.json").write_text(json.dumps(manifest))
    ck.write_cohort_checksums(cs)   # real blake2b over the artefacts above
    return cs


def _fake_reader(chunk_dir):
    return np.array([[20.0, 8.0], [8.0, 10.0]])


def test_decrypted_cohort_accumulates_with_checksum_verification(tmp_path):
    wrapper = _make_gpg_wrapper(tmp_path)
    env = _env(wrapper)
    cohort = _build_real_cohort(tmp_path, study="cohortA")
    upload = tmp_path / "upload"
    rebuilt = tmp_path / "rebuilt"
    assert _encrypt(env, cohort, upload, study="cohortA").returncode == 0
    assert _decrypt(env, upload, rebuilt, study="cohortA").returncode == 0

    # checksums.json round-tripped, so accumulate verifies it against the rebuilt tree
    precursor = tmp_path / "precursor"
    agg.accumulate(rebuilt, precursor, chunk_reader=_fake_reader,
                   pair_batch_rows=16, verify_checksums=True)
    pm = agg.read_precursor_manifest(precursor)
    assert pm["n_cohorts"] == 1
    assert pm["cohorts"][0]["study_name"] == "cohortA"
```

- [ ] **Step 2: Run the test to verify it passes**

Run: `pytest tests/genetics/test_ld_decrypt_cohort.py::test_decrypted_cohort_accumulates_with_checksum_verification -v`
Expected: PASS — `accumulate` re-verifies the round-tripped `checksums.json` against the rebuilt tree and folds the cohort in. (If this fails with a `ChecksumError`, the round-trip is not byte-exact — investigate before proceeding; do not disable verification.)

- [ ] **Step 3: Run the full decrypt + aggregate suites as a regression check**

Run: `pytest tests/genetics/test_ld_decrypt_cohort.py tests/genetics/test_ld_aggregate.py tests/genetics/test_ld_encrypt_cohort.py -v`
Expected: PASS — all decrypt, aggregate, and encrypt tests.

- [ ] **Step 4: Commit**

```bash
git add tests/genetics/test_ld_decrypt_cohort.py
git commit -m "test(15): end-to-end decrypt -> accumulate with checksum verification"
```

---

## Notes for the implementer

- **Why manifest.json moves last in the scaffold restore:** the resume guard treats `output_dir/manifest.json` as the "scaffold complete" marker. Moving it last means an interrupted scaffold restore is correctly re-done on the next run.
- **Why chunks publish via `mv` of the final dir:** `tar` extracts into a temp dir, then the complete `chr<C>/chunk_<N>` directory is renamed into place in one step, so the chunk's resume marker (`-d "${target}"`) only ever sees complete chunks.
- **md5sum semantics:** the `.md5sum` files store a bare filename (`<base>.tgz`), so `md5sum -c` must run with cwd = the dir holding the decrypted tar. That is why `input_dir`/`output_dir` are absolutised up front and verification runs inside `${staging}`.
- **Real gpg in tests:** the wrapper uses an isolated `--homedir` and a loopback passphrase, so tests need a working `gpg` binary but never touch the user's keyring. This matches `test_ld_encrypt_cohort.py`.
- **Operator workflow (documented in the spec):** decrypt cohort A → `bash 15b-ld_aggregate_stats.sh` (or `ld_aggregate_stats.py --mode accumulate --cohort-dir <A>`) → repeat per cohort → finalise. `${ld_prepare_dir}` is the single `--cohort-dir`; point it at each decrypt `output_dir` in turn.
```
