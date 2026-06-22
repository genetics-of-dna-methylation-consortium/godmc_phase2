# Section 15 (LD) Cohort Run-Upload-Delete Wrapper Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give a cohort one resumable command (`15c-ld_run_upload.sh`) that produces, ships to Imperial, and locally reclaims their section-15 LD outputs one chromosome at a time, so they need ~150 GB free disk instead of ~2 TB.

**Architecture:** A new bash orchestrator drives the existing per-chromosome compute worker `ld_prepare_stats.py` once per chromosome into its own output directory, then streams each `A_blocks` chunk through encrypt → rsync → verify → delete before moving on. Deletion is gated on a verified upload; a per-chromosome `.uploaded` sentinel makes the whole loop resumable. The wrapper exposes override hooks (`LD_PREPARE_CMD`, `LD_SHIP_CMD`, `GPG`) so every behaviour is testable without Hail, gpg-agent, or the network.

**Tech Stack:** bash, GnuPG (symmetric AES256), `tar`/`gzip`/`md5sum`, rsync-over-SSH, pytest + `subprocess`.

**Spec:** `docs/superpowers/specs/2026-06-22-section15-cohort-run-upload-wrapper-design.md`

## Global Constraints

- Encryption command, verbatim (matches all other sections): `gpg --output X.aes --symmetric --cipher-algo AES256 X`. The gpg binary is read from `${GPG:-gpg}` so tests can inject a non-interactive wrapper.
- `md5sum` is computed on the **plaintext** `.tgz` before encryption; the `.tgz` is removed immediately after encrypting to bound peak disk.
- Per-chromosome compute writes to `results/15/cohort_stats/chr<C>/` (its own directory). `ld_prepare_stats.py` is NOT modified — it already accepts `--chromosome` and `--output-dir`.
- Encrypted artifacts stage in `${section_15_dir}/upload/` = `results/15/upload/`.
- Artifact naming is chromosome-namespaced via a per-chr study tag `${study_name}_chr<C>`: scaffold `${study_name}_chr<C>_15_scaffold.tgz.aes` (+`.md5sum`); per-chunk `${study_name}_chr<C>_15_chr<C>_chunk_<N>.tgz.aes` (+`.md5sum`).
- **Nothing is deleted until its upload is verified.** A chromosome's `A_blocks` is removed only after every chunk and the scaffold have shipped successfully and the `.uploaded` sentinel is written.
- Transport is rsync-over-SSH with `--partial --append --checksum`, isolated in one function (`ship_file`) and overridable via `LD_SHIP_CMD` so the concrete Imperial endpoint can be swapped later. Endpoint config: fixed `imperial_host`/`imperial_path` in `resources/parameters`; per-cohort `imperial_user`/`imperial_key` in `config`.
- Cohort-specific credentials go in `config` (like `sftp_username`/`key`); fixed consortium endpoints go in `resources/parameters` (like `sftp_address`/`sftp_path`).
- The wrapper does NOT use `check_upload.sh` (that targets Bristol SFTP); section 15 ships to Imperial.
- Central reassembly/decryption is OUT OF SCOPE (separate follow-up spec).

---

### Task 1: Config + parameters — Imperial endpoint vars and `ld_chromosomes`

**Files:**
- Modify: `config.example` (add `imperial_user`, `imperial_key` near the existing `sftp_username`/`key` block)
- Modify: `resources/parameters` (add `imperial_host`, `imperial_path`, `ld_chromosomes` near the existing `ld_*` block ~line 262 and the sftp block ~line 285)
- Test: `tests/genetics/test_ld_run_upload_config.py`

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces: config vars `imperial_user`, `imperial_key` (cohort-set, default empty); `imperial_host`, `imperial_path` (fixed, placeholder until Imperial finalised); `ld_chromosomes` (default the autosome list `1 2 3 ... 22`, space-separated, overridable via env).

- [ ] **Step 1: Write the failing test**

Create `tests/genetics/test_ld_run_upload_config.py`:

```python
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_config_example_has_imperial_creds():
    text = (REPO_ROOT / "config.example").read_text()
    assert "imperial_user=" in text
    assert "imperial_key=" in text


def test_parameters_has_imperial_endpoint_and_chromosomes():
    text = (REPO_ROOT / "resources" / "parameters").read_text()
    assert "imperial_host=" in text
    assert "imperial_path=" in text
    # default autosome list, env-overridable
    assert 'ld_chromosomes="${ld_chromosomes:-' in text
    # first and last autosome present in the default
    assert "1 2 3" in text and "22" in text
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload_config.py -v`
Expected: FAIL — the new vars are not present yet.

- [ ] **Step 3: Add the cohort credential vars to `config.example`**

In `config.example`, immediately after the `key=""` line (the existing SFTP key, ~line 18), insert:

```bash

# Imperial LD-upload (section 15) account, provided by the GoDMC developers.
imperial_user=""
# Path to the SSH private key for the Imperial LD-upload account (often the same as `key`).
imperial_key=""
```

- [ ] **Step 4: Add the fixed endpoint + chromosome list to `resources/parameters`**

In `resources/parameters`, immediately after the `ld_min_cohorts=...` line (~line 262), insert:

```bash
ld_chromosomes="${ld_chromosomes:-1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22}"
```

Then, in the SFTP details block (after the `sftp_address=...`/`sftp_path` lines, ~line 309), insert:

```bash

# Imperial section-15 (LD) rsync endpoint. Placeholder until the endpoint is
# finalised; the wrapper's transport step is isolated so this can change in one place.
imperial_host="${imperial_host:-TODO.imperial.ac.uk}"
imperial_path="${imperial_path:-/godmc/ld_upload}"
```

> Note: `TODO.imperial.ac.uk` is the agreed placeholder value for an undecided endpoint (it is a default, not a plan placeholder); the wrapper never contacts it in tests because `LD_SHIP_CMD` overrides the transport.

- [ ] **Step 5: Run test to verify it passes**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload_config.py -v`
Expected: PASS (2 tests).

- [ ] **Step 6: Commit**

```bash
git add config.example resources/parameters tests/genetics/test_ld_run_upload_config.py
git commit -m "feat(15): Imperial LD-upload endpoint vars + ld_chromosomes default"
```

---

### Task 2: Wrapper scaffolding — resume sentinel + transport indirection

**Files:**
- Create: `15c-ld_run_upload.sh`
- Test: `tests/genetics/test_ld_run_upload.py`

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces (functions, available when the script is sourced):
  - `chr_done <outdir>` — exit 0 iff `<outdir>/.uploaded` exists.
  - `ship_file <path>` — ships one file; if `${LD_SHIP_CMD}` is set, runs `"${LD_SHIP_CMD}" <path>` and returns its exit code; otherwise rsyncs to Imperial. Returns the transport's exit code.
- The script must define functions at top level and only run `main` when executed directly (so tests can source it without side effects): guard with `if [ "${BASH_SOURCE[0]}" = "${0}" ]; then main "$@"; fi`.

- [ ] **Step 1: Write the failing test**

Create `tests/genetics/test_ld_run_upload.py`:

```python
import os
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "15c-ld_run_upload.sh"


def _call(func_call, env=None, extra=""):
    """Source the wrapper (main does not run) and invoke one shell snippet."""
    full_env = dict(os.environ)
    if env:
        full_env.update(env)
    return subprocess.run(
        ["bash", "-c", f'source "{SCRIPT}"; {extra} {func_call}'],
        capture_output=True, text=True, env=full_env,
    )


def test_parses():
    r = subprocess.run(["bash", "-n", str(SCRIPT)], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr


def test_chr_done_detects_sentinel(tmp_path):
    outdir = tmp_path / "chr22"
    outdir.mkdir()
    # no sentinel yet
    assert _call(f'chr_done "{outdir}"').returncode != 0
    (outdir / ".uploaded").touch()
    assert _call(f'chr_done "{outdir}"').returncode == 0


def test_ship_file_uses_override_and_propagates_exit(tmp_path):
    f = tmp_path / "artifact.aes"
    f.write_text("payload")
    dest = tmp_path / "remote"
    dest.mkdir()
    ok = tmp_path / "ship_ok.sh"
    ok.write_text(f'#!/usr/bin/env bash\ncp "$1" "{dest}/"\n')
    ok.chmod(0o755)
    r = _call(f'ship_file "{f}"', env={"LD_SHIP_CMD": str(ok)})
    assert r.returncode == 0, r.stderr
    assert (dest / "artifact.aes").read_text() == "payload"

    fail = tmp_path / "ship_fail.sh"
    fail.write_text("#!/usr/bin/env bash\nexit 7\n")
    fail.chmod(0o755)
    assert _call(f'ship_file "{f}"', env={"LD_SHIP_CMD": str(fail)}).returncode == 7
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload.py -v`
Expected: FAIL — `15c-ld_run_upload.sh` does not exist.

- [ ] **Step 3: Create the wrapper with scaffolding + the two functions**

Create `15c-ld_run_upload.sh`:

```bash
#!/usr/bin/env bash
# 15c-ld_run_upload.sh — cohort-side section-15 (LD) run -> encrypt -> upload
# (Imperial) -> verify -> delete, one chromosome at a time, so peak disk stays
# at roughly one chromosome's A_blocks (~100 GB) instead of the ~1.3 TB a
# whole-genome run would require.
#
# Usage: ./15c-ld_run_upload.sh -c <config>
#
# Resumable: each chromosome writes results/15/cohort_stats/chr<C>/.uploaded
# only after every artifact has shipped; re-running skips completed chromosomes.
#
# Test/override hooks (no effect in normal operation):
#   GPG            gpg binary/wrapper (default: gpg)
#   LD_SHIP_CMD    command taking one file path; replaces the rsync transport
#   LD_PREPARE_CMD command "<chr> <outdir>"; replaces the 15a compute step
set -euo pipefail

GPG="${GPG:-gpg}"

# chr_done <outdir> — true if the chromosome's upload sentinel exists.
chr_done () {
	[ -f "$1/.uploaded" ]
}

# ship_file <path> — transfer one file to Imperial (or via the test override).
ship_file () {
	local path="$1"
	if [ -n "${LD_SHIP_CMD:-}" ]; then
		"${LD_SHIP_CMD}" "${path}"
		return $?
	fi
	rsync --partial --append --checksum \
		-e "ssh -i ${imperial_key}" \
		"${path}" "${imperial_user}@${imperial_host}:${imperial_path}/"
}

main () {
	source resources/setup.sh "$@"
	set -- $concatenated
	echo "[15c] main is implemented in a later task"
}

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
	main "$@"
fi
```

Make it executable:

```bash
chmod +x 15c-ld_run_upload.sh
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload.py -v`
Expected: PASS (3 tests).

- [ ] **Step 5: Commit**

```bash
git add 15c-ld_run_upload.sh tests/genetics/test_ld_run_upload.py
git commit -m "feat(15): 15c wrapper scaffold — resume sentinel + transport indirection"
```

---

### Task 3: Per-chunk encrypt → ship → delete-on-success

**Files:**
- Modify: `15c-ld_run_upload.sh`
- Test: `tests/genetics/test_ld_run_upload.py`

**Interfaces:**
- Consumes: Task 2's `ship_file`, `${GPG}`.
- Produces:
  - `encrypt_archive <src_parent> <member> <out_dir> <base>` — `tar czf out_dir/base.tgz -C src_parent member`, `md5sum` → `out_dir/base.md5sum`, `gpg` → `out_dir/base.tgz.aes`, then `rm` the `.tgz`.
  - `process_chunk <chunk_dir> <ablocks_dir> <out_dir> <study_tag>` — encrypts the chunk (base `${study_tag}_15_<chr>_<chunk>`), ships the `.aes` and `.md5sum`; on success removes the source chunk dir AND the local `.aes`/`.md5sum` and returns 0; on ship failure leaves everything in place and returns 1.

- [ ] **Step 1: Write the failing test**

Append to `tests/genetics/test_ld_run_upload.py`:

```python
PASSPHRASE = "testpass"


def _gpg_wrapper(tmp_path):
    home = tmp_path / "gnupg"
    home.mkdir(mode=0o700)
    w = tmp_path / "gpg_batch.sh"
    w.write_text(
        "#!/usr/bin/env bash\n"
        f'exec gpg --homedir "{home}" --batch --yes '
        f'--pinentry-mode loopback --passphrase "{PASSPHRASE}" "$@"\n'
    )
    w.chmod(0o755)
    return w


def _ship_ok(tmp_path):
    dest = tmp_path / "remote"
    dest.mkdir(exist_ok=True)
    s = tmp_path / "ship_ok.sh"
    s.write_text(f'#!/usr/bin/env bash\ncp "$1" "{dest}/"\n')
    s.chmod(0o755)
    return s, dest


def _make_chunk(ablocks, chrom="chr22", chunk="chunk_0"):
    d = ablocks / chrom / chunk
    d.mkdir(parents=True)
    (d / "part-00000").write_bytes(b"blockmatrix-bytes")
    (d / "metadata.json").write_text('{"block":1}')
    return d


def test_process_chunk_success_ships_and_deletes(tmp_path):
    ablocks = tmp_path / "A_blocks"
    chunk = _make_chunk(ablocks)
    out = tmp_path / "upload"; out.mkdir()
    gpg = _gpg_wrapper(tmp_path)
    ship, dest = _ship_ok(tmp_path)
    r = _call(
        f'process_chunk "{chunk}" "{ablocks}" "{out}" "study_chr22"',
        env={"GPG": str(gpg), "LD_SHIP_CMD": str(ship)},
    )
    assert r.returncode == 0, r.stderr
    base = "study_chr22_15_chr22_chunk_0"
    # shipped to the fake remote
    assert (dest / f"{base}.tgz.aes").is_file()
    assert (dest / f"{base}.md5sum").is_file()
    # local source + staged artifacts removed
    assert not chunk.exists()
    assert list(out.glob("*")) == []


def test_process_chunk_failed_ship_keeps_everything(tmp_path):
    ablocks = tmp_path / "A_blocks"
    chunk = _make_chunk(ablocks)
    out = tmp_path / "upload"; out.mkdir()
    gpg = _gpg_wrapper(tmp_path)
    fail = tmp_path / "ship_fail.sh"
    fail.write_text("#!/usr/bin/env bash\nexit 1\n")
    fail.chmod(0o755)
    r = _call(
        f'process_chunk "{chunk}" "{ablocks}" "{out}" "study_chr22"',
        env={"GPG": str(gpg), "LD_SHIP_CMD": str(fail)},
    )
    assert r.returncode == 1
    # source chunk NOT deleted on a failed upload
    assert chunk.exists()
    assert (chunk / "part-00000").is_file()
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload.py -k "process_chunk" -v`
Expected: FAIL — `process_chunk`/`encrypt_archive` are not defined.

- [ ] **Step 3: Add the two functions**

In `15c-ld_run_upload.sh`, insert these functions immediately after `ship_file` (before `main`):

```bash
# encrypt_archive <src_parent> <member> <out_dir> <base>
# tar+md5(plaintext)+gpg one path; removes the intermediate .tgz.
encrypt_archive () {
	local src_parent="$1" member="$2" out_dir="$3" base="$4"
	mkdir -p "${out_dir}"
	tar czf "${out_dir}/${base}.tgz" -C "${src_parent}" "${member}"
	( cd "${out_dir}" && md5sum "${base}.tgz" > "${base}.md5sum" )
	"${GPG}" --output "${out_dir}/${base}.tgz.aes" \
		--symmetric --cipher-algo AES256 "${out_dir}/${base}.tgz"
	rm -f "${out_dir}/${base}.tgz"
}

# process_chunk <chunk_dir> <ablocks_dir> <out_dir> <study_tag>
# encrypt -> ship .aes + .md5sum -> on success delete source + staged artifacts;
# on ship failure leave everything in place and return 1.
process_chunk () {
	local chunk_dir="$1" ablocks_dir="$2" out_dir="$3" study_tag="$4"
	local chr_name chunk_name base
	chr_name="$(basename "$(dirname "${chunk_dir}")")"
	chunk_name="$(basename "${chunk_dir}")"
	base="${study_tag}_15_${chr_name}_${chunk_name}"
	encrypt_archive "${ablocks_dir}" "${chr_name}/${chunk_name}" "${out_dir}" "${base}"
	if ship_file "${out_dir}/${base}.tgz.aes" && ship_file "${out_dir}/${base}.md5sum"; then
		rm -rf "${chunk_dir}" "${out_dir}/${base}.tgz.aes" "${out_dir}/${base}.md5sum"
		echo "[15c] shipped ${base}"
		return 0
	fi
	echo "[15c] ship FAILED for ${base}; leaving source in place" >&2
	return 1
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload.py -v`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add 15c-ld_run_upload.sh tests/genetics/test_ld_run_upload.py
git commit -m "feat(15): per-chunk encrypt+ship with delete gated on verified upload"
```

---

### Task 4: Per-chromosome processing + main loop

**Files:**
- Modify: `15c-ld_run_upload.sh`
- Test: `tests/genetics/test_ld_run_upload.py`

**Interfaces:**
- Consumes: Task 2's `chr_done`/`ship_file`, Task 3's `encrypt_archive`/`process_chunk`, `${GPG}`.
- Produces:
  - `prepare_chromosome <C> <outdir>` — runs the 15a compute for chromosome `C` into `<outdir>`; if `${LD_PREPARE_CMD}` is set, runs `"${LD_PREPARE_CMD}" <C> <outdir>` instead (test/override).
  - `check_chromosome <outdir>` — returns 0 iff `manifest.json`, `variants.tsv.gz`, `D.npy`, `B.npy` and the `A_blocks/` dir are present; returns 1 otherwise.
  - `ship_scaffold <outdir> <out_dir> <study_tag>` — encrypts the small scaffold files into `${study_tag}_15_scaffold`, ships `.aes` + `.md5sum`; returns 0 on success (and removes the staged artifacts), 1 on failure.
  - `process_chromosome <C> <study> <cohort_root> <out_dir>` — skip if done; prepare; check; ship every chunk; ship scaffold; on full success write `<cohort_root>/chr<C>/.uploaded` and `rm -rf` that chromosome's `A_blocks`. Returns non-zero if any step fails (leaving data intact).
  - `main` — sources config, iterates `${ld_chromosomes}`, calls `process_chromosome` per chr, continues past a failed chromosome, and exits non-zero if any chromosome failed.

- [ ] **Step 1: Write the failing test**

Append to `tests/genetics/test_ld_run_upload.py`:

```python
def _prepare_stub(tmp_path):
    """A fake 15a: creates scaffold files + two A_blocks chunks for <chr> in <outdir>."""
    s = tmp_path / "prepare_stub.sh"
    s.write_text(
        "#!/usr/bin/env bash\n"
        'set -euo pipefail\n'
        'C="$1"; OUT="$2"\n'
        'mkdir -p "$OUT/A_blocks/chr${C}/chunk_0" "$OUT/A_blocks/chr${C}/chunk_1"\n'
        'echo data0 > "$OUT/A_blocks/chr${C}/chunk_0/part-00000"\n'
        'echo data1 > "$OUT/A_blocks/chr${C}/chunk_1/part-00000"\n'
        'echo "{}" > "$OUT/manifest.json"\n'
        'printf "" | gzip > "$OUT/variants.tsv.gz"\n'
        'echo D > "$OUT/D.npy"; echo B > "$OUT/B.npy"\n'
        'echo "{}" > "$OUT/checksums.json"; echo qc > "$OUT/qc_report.txt"\n'
    )
    s.chmod(0o755)
    return s


def _process_chr_env(tmp_path):
    gpg = _gpg_wrapper(tmp_path)
    ship, dest = _ship_ok(tmp_path)
    prep = _prepare_stub(tmp_path)
    return {
        "GPG": str(gpg), "LD_SHIP_CMD": str(ship), "LD_PREPARE_CMD": str(prep),
    }, dest


def test_process_chromosome_success_deletes_ablocks_and_marks(tmp_path):
    root = tmp_path / "cohort_stats"; root.mkdir()
    out = tmp_path / "upload"; out.mkdir()
    env, dest = _process_chr_env(tmp_path)
    r = _call(f'process_chromosome 22 "study" "{root}" "{out}"', env=env)
    assert r.returncode == 0, r.stderr
    chrdir = root / "chr22"
    assert (chrdir / ".uploaded").is_file()          # marked
    assert not (chrdir / "A_blocks").exists()         # big data reclaimed
    assert (chrdir / "D.npy").is_file()               # scaffold kept locally
    # scaffold + 2 chunks shipped (each: .tgz.aes + .md5sum)
    assert (dest / "study_chr22_15_scaffold.tgz.aes").is_file()
    assert (dest / "study_chr22_15_chr22_chunk_0.tgz.aes").is_file()
    assert (dest / "study_chr22_15_chr22_chunk_1.tgz.aes").is_file()


def test_process_chromosome_skips_when_already_done(tmp_path):
    root = tmp_path / "cohort_stats"; (root / "chr22").mkdir(parents=True)
    (root / "chr22" / ".uploaded").touch()
    out = tmp_path / "upload"; out.mkdir()
    env, dest = _process_chr_env(tmp_path)
    r = _call(f'process_chromosome 22 "study" "{root}" "{out}"', env=env)
    assert r.returncode == 0, r.stderr
    assert "skip" in r.stdout.lower()
    assert list(dest.glob("*")) == []   # nothing shipped


def test_process_chromosome_failed_upload_keeps_ablocks_no_sentinel(tmp_path):
    root = tmp_path / "cohort_stats"; root.mkdir()
    out = tmp_path / "upload"; out.mkdir()
    gpg = _gpg_wrapper(tmp_path)
    prep = _prepare_stub(tmp_path)
    fail = tmp_path / "ship_fail.sh"
    fail.write_text("#!/usr/bin/env bash\nexit 1\n")
    fail.chmod(0o755)
    env = {"GPG": str(gpg), "LD_SHIP_CMD": str(fail), "LD_PREPARE_CMD": str(prep)}
    r = _call(f'process_chromosome 22 "study" "{root}" "{out}"', env=env)
    assert r.returncode != 0
    chrdir = root / "chr22"
    assert not (chrdir / ".uploaded").exists()        # not marked
    assert (chrdir / "A_blocks").exists()             # data preserved
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload.py -k "process_chromosome" -v`
Expected: FAIL — `process_chromosome` and helpers are not defined.

- [ ] **Step 3: Add the per-chromosome functions and the real `main`**

In `15c-ld_run_upload.sh`, insert these functions immediately after `process_chunk` (before `main`):

```bash
# prepare_chromosome <C> <outdir> — run 15a compute for chromosome C into outdir.
prepare_chromosome () {
	local chr="$1" outdir="$2"
	mkdir -p "${outdir}"
	if [ -n "${LD_PREPARE_CMD:-}" ]; then
		"${LD_PREPARE_CMD}" "${chr}" "${outdir}"
		return $?
	fi
	local hail_runtime_args=""
	if [ -n "${ld_hail_local_cores:-}" ]; then
		hail_runtime_args="${hail_runtime_args} --hail-local-cores ${ld_hail_local_cores}"
	fi
	if [ -n "${ld_hail_driver_memory_gb:-}" ]; then
		hail_runtime_args="${hail_runtime_args} --hail-driver-memory-gb ${ld_hail_driver_memory_gb}"
	fi
	python "${scripts_directory}/resources/genetics/ld_prepare_stats.py" \
		--study-name "${study_name}" \
		--bfile "${bfile}" \
		--covariates "${covariates_intersect}" \
		--output-dir "${outdir}" \
		--log-file "${section_15a_logfile}" \
		--chromosome "${chr}" \
		--hail-partitions "${ld_hail_partitions}" \
		--hail-tmp-dir "${ld_hail_tmp_dir}" \
		--a-block-size "${ld_a_block_size}" \
		--a-chunk-rows "${ld_a_chunk_rows}" \
		--a-max-dense-gb "${ld_a_max_dense_gb}" \
		${hail_runtime_args}
}

# check_chromosome <outdir> — required scaffold + A_blocks present.
check_chromosome () {
	local outdir="$1" f
	for f in manifest.json variants.tsv.gz D.npy B.npy; do
		if [ ! -f "${outdir}/${f}" ]; then
			echo "[15c] check failed: missing ${f} in ${outdir}" >&2
			return 1
		fi
	done
	if [ ! -d "${outdir}/A_blocks" ]; then
		echo "[15c] check failed: missing A_blocks in ${outdir}" >&2
		return 1
	fi
}

# ship_scaffold <outdir> <out_dir> <study_tag>
ship_scaffold () {
	local outdir="$1" out_dir="$2" study_tag="$3"
	local base="${study_tag}_15_scaffold"
	mkdir -p "${out_dir}"
	tar czf "${out_dir}/${base}.tgz" -C "${outdir}" \
		manifest.json variants.tsv.gz D.npy B.npy checksums.json qc_report.txt
	( cd "${out_dir}" && md5sum "${base}.tgz" > "${base}.md5sum" )
	"${GPG}" --output "${out_dir}/${base}.tgz.aes" \
		--symmetric --cipher-algo AES256 "${out_dir}/${base}.tgz"
	rm -f "${out_dir}/${base}.tgz"
	if ship_file "${out_dir}/${base}.tgz.aes" && ship_file "${out_dir}/${base}.md5sum"; then
		rm -f "${out_dir}/${base}.tgz.aes" "${out_dir}/${base}.md5sum"
		echo "[15c] shipped ${base}"
		return 0
	fi
	echo "[15c] ship FAILED for ${base}; leaving source in place" >&2
	return 1
}

# process_chromosome <C> <study> <cohort_root> <out_dir>
process_chromosome () {
	local chr="$1" study="$2" cohort_root="$3" out_dir="$4"
	local outdir="${cohort_root}/chr${chr}"
	if chr_done "${outdir}"; then
		echo "[15c] chr${chr} already uploaded; skip"
		return 0
	fi
	prepare_chromosome "${chr}" "${outdir}"
	check_chromosome "${outdir}" || return 1
	local ablocks="${outdir}/A_blocks" tag="${study}_chr${chr}"
	local chunk_dir
	shopt -s nullglob
	for chunk_dir in "${ablocks}"/chr*/chunk_*; do
		[ -d "${chunk_dir}" ] || continue
		if ! process_chunk "${chunk_dir}" "${ablocks}" "${out_dir}" "${tag}"; then
			shopt -u nullglob
			return 1
		fi
	done
	shopt -u nullglob
	ship_scaffold "${outdir}" "${out_dir}" "${tag}" || return 1
	touch "${outdir}/.uploaded"
	rm -rf "${ablocks}"
	echo "[15c] chr${chr} complete"
}
```

Then REPLACE the placeholder `main` body with the real loop. Replace:

```bash
main () {
	source resources/setup.sh "$@"
	set -- $concatenated
	echo "[15c] main is implemented in a later task"
}
```

with:

```bash
main () {
	source resources/setup.sh "$@"
	set -- $concatenated

	exec &> >(tee "${section_15a_logfile}")
	print_version

	local cohort_root="${ld_prepare_dir}"
	local upload_dir="${section_15_dir}/upload"
	mkdir -p "${cohort_root}" "${upload_dir}" "${ld_hail_tmp_dir}"

	if [ ! -f "${bfile}.bed" ]; then
		echo "Problem: cleaned section-02 genotype files are required at ${bfile}"
		exit 1
	fi
	if [ ! -f "${covariates_intersect}" ]; then
		echo "Problem: section-03a mQTL-aligned covariates are required at ${covariates_intersect}"
		exit 1
	fi

	local failed=0 chr
	for chr in ${ld_chromosomes}; do
		echo "[15c] ===== chromosome ${chr} ====="
		if ! process_chromosome "${chr}" "${study_name}" "${cohort_root}" "${upload_dir}"; then
			echo "[15c] chromosome ${chr} did NOT complete; re-run to resume" >&2
			failed=1
		fi
	done

	if [ "${failed}" -ne 0 ]; then
		echo "[15c] one or more chromosomes failed; re-run ./15c-ld_run_upload.sh to resume"
		exit 1
	fi
	echo "Successfully ran and uploaded all section-15 LD cohort chromosomes"
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload.py -v`
Expected: PASS (8 tests).

- [ ] **Step 5: Commit**

```bash
git add 15c-ld_run_upload.sh tests/genetics/test_ld_run_upload.py
git commit -m "feat(15): per-chromosome run->ship->reclaim loop + resumable main"
```

---

### Task 5: Cohort runbook (wiki)

**Files:**
- Rewrite: `godmc_phase2.wiki/Run-federated-LD-reference-panel.md`
- Test: `tests/genetics/test_ld_wiki_runbook.py`

**Interfaces:**
- Consumes: the finished `15c-ld_run_upload.sh` and the config vars from Task 1.
- Produces: a cohort operator runbook. No code; the "test" asserts the page documents the real command and no longer carries the stale "not yet implemented" status.

- [ ] **Step 1: Write the failing test**

Create `tests/genetics/test_ld_wiki_runbook.py`:

```python
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
WIKI = REPO_ROOT / "godmc_phase2.wiki" / "Run-federated-LD-reference-panel.md"


def test_runbook_documents_wrapper_and_is_current():
    text = WIKI.read_text()
    # documents the one cohort command
    assert "15c-ld_run_upload.sh" in text
    # documents the disk expectation that motivates the design
    assert "150 GB" in text
    # documents the Imperial endpoint config the operator must set
    assert "imperial_user" in text and "imperial_key" in text
    # documents resumability
    assert ".uploaded" in text
    # stale status removed
    assert "not yet implemented" not in text
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_wiki_runbook.py -v`
Expected: FAIL — the page still says "not yet implemented" and does not mention the wrapper.

- [ ] **Step 3: Rewrite the wiki page**

Replace the entire contents of `godmc_phase2.wiki/Run-federated-LD-reference-panel.md` with:

```markdown
# MODULE STATUS

**Developers:** Toby Clark

**Scripts status:** Cohort-side complete (run + upload). Central aggregation (15b) in development.

**Prerequisite scripts:** `00-setup_folders.sh`, `01-check_data.sh`, `02a-snp_data.sh`, `03a-methylation_variables.sh`

**Data upload method:** `15c-ld_run_upload.sh` — per-chromosome run → encrypt → upload (Imperial) → delete.

## Overview

Section 15 builds a federated LD reference panel from cohort-level summary
statistics (`A`, `B`, `D` matrices) rather than centrally stored individual
genotypes. The first release is autosomes only and uses the cleaned section-02
genotypes at `${bfile}`, aligned to the mQTL sample set from
`03a-methylation_variables.sh`.

Cohorts run **one command**, `15c-ld_run_upload.sh`. Central aggregation
(`15b-ld_aggregate_stats.sh`) is run by the developers, not by cohorts.

## Why section 15 has its own run+upload command

Every other module produces small results, so the convention is "run the
script, then `check_upload.sh <n> check` / `<n> upload`". Section 15 cannot do
that: its genome-wide `A_blocks` output is ~1.3 TB per cohort, and that disk is
consumed *while computing*, before any upload could start. So section 15
processes **one chromosome at a time**: it computes a chromosome, encrypts and
uploads it, verifies the transfer, then deletes that chromosome's large
`A_blocks` before moving to the next. **Peak local disk is roughly one
chromosome (~100 GB), so you need about 150 GB free — not 2 TB.**

## Prerequisites

- Completed: `00-setup_folders.sh`, `01-check_data.sh`, `02a-snp_data.sh`,
  `03a-methylation_variables.sh`.
- The Hail environment from `resources/genetics/hail_env.yml`.
- **~150 GB free disk** under your results directory.
- Imperial upload credentials (provided by the GoDMC developers) and the LD
  encryption passphrase (shared out-of-band).

## Configuration

In your `config` file set the Imperial upload account:

    imperial_user="<provided by developers>"
    imperial_key="~/.ssh/id_rsa"   # SSH key for the Imperial account

Optional overrides (environment variables; defaults are fine for production):

- `ld_chromosomes` — chromosomes to process, default `1 2 … 22`. Set to a single
  value (e.g. `22`) for a pilot.
- `ld_hail_local_cores`, `ld_hail_driver_memory_gb` — Hail/Spark resources.

## Running

    ./15c-ld_run_upload.sh -c <config>

For each chromosome the script: computes the LD statistics into
`results/15/cohort_stats/chr<C>/`, checks the outputs, encrypts each `A_blocks`
chunk and the small scaffold (`manifest.json`, `variants.tsv.gz`, `D.npy`,
`B.npy`) with symmetric AES256, uploads them to Imperial, verifies the transfer,
then deletes that chromosome's `A_blocks`. You are prompted once for the
encryption passphrase (gpg-agent caches it for the session).

## Resumability and safety

- **Safe to re-run.** Each chromosome writes `results/15/cohort_stats/chr<C>/.uploaded`
  only after everything for it has shipped. Re-running skips completed
  chromosomes, so an interrupted run (lost connection, time limit) resumes where
  it stopped — just run the same command again.
- **Nothing is deleted until its upload is verified.** A chromosome's `A_blocks`
  is removed only after every chunk and the scaffold have transferred
  successfully. An interrupted chromosome keeps all its data and is recomputed
  on the next run.
- **Progress:** the chromosomes with a `.uploaded` file under
  `results/15/cohort_stats/` are the ones already shipped.

## Expectations

On the 1000G pilot (8 cores / 32 GB), chr1 took ≈ 2.5 h and peaked at ~100 GB of
`A_blocks`; genome-wide is ≈ 33 h streamed. Your hardware will differ — treat
these as a reference, not a guarantee.

## Covariate schema

The first release is **intercept only**. The `03a` intersected covariate file is
used only to align the LD sample set to the mQTL dataset; its covariate columns
(age, sex, PCs, cell counts, smoking, batch, methylation PCs) are NOT
residualised from genotypes — they are phenotype-side adjustments in the mQTL
pipeline. Any future non-intercept schema must be an explicitly versioned panel
specification matching the SNP-side association model.
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_wiki_runbook.py -v`
Expected: PASS (1 test).

- [ ] **Step 5: Full section-15 suite green + commit**

Run: `cd /mnt/data2/tobyc/godmc_phase2 && python -m pytest tests/genetics/test_ld_run_upload.py tests/genetics/test_ld_run_upload_config.py tests/genetics/test_ld_wiki_runbook.py -v`
Expected: PASS (11 tests).

```bash
git add godmc_phase2.wiki/Run-federated-LD-reference-panel.md tests/genetics/test_ld_wiki_runbook.py
git commit -m "docs(15): cohort runbook for 15c per-chromosome run-upload-delete"
```

> The `godmc_phase2.wiki/` directory is its own git repo (GitHub wiki). The commit above stages the file in the main repo's working tree; to publish the wiki, also commit + push inside `godmc_phase2.wiki/` (`cd godmc_phase2.wiki && git add … && git commit && git push`). Confirm with the user before pushing the wiki.

---

## Self-Review

**1. Spec coverage:**
- New `15c-ld_run_upload.sh` orchestrator, one command → Tasks 2–4. ✓
- Per-chromosome loop, skip-if-done sentinel, resumable → Task 2 (`chr_done`), Task 4 (`main`, `process_chromosome`). ✓
- Compute via existing `ld_prepare_stats.py` per-chr `--output-dir` (no worker change) → Task 4 (`prepare_chromosome`). ✓
- Streamed per-chunk encrypt → upload → verify → delete, disk bounded to ~one chromosome → Task 3 (`process_chunk`). ✓
- Delete only after verified upload → Tasks 3–4 (ship-gated `rm`, sentinel after success). ✓
- Imperial rsync transport, isolated + swappable (`LD_SHIP_CMD`), new config vars → Task 1, Task 2 (`ship_file`). ✓
- Chromosome-namespaced artifact names (`${study_name}_chr<C>`), no scaffold collision → Tasks 3–4. ✓
- Encryption convention (`gpg --symmetric --cipher-algo AES256`, plaintext md5, rm .tgz) → Task 3 (`encrypt_archive`), Task 4 (`ship_scaffold`). ✓
- Wiki rewritten as cohort runbook; stale headers fixed → Task 5. ✓
- Central reassembly explicitly out of scope → not in plan (correct). ✓

**2. Placeholder scan:** No TBD/TODO as plan instructions. The only literal "TODO" is the agreed `imperial_host` default value (`TODO.imperial.ac.uk`), explicitly flagged as an intentional placeholder default, never contacted in tests. ✓

**3. Type/name consistency:** Function names (`chr_done`, `ship_file`, `encrypt_archive`, `process_chunk`, `prepare_chromosome`, `check_chromosome`, `ship_scaffold`, `process_chromosome`, `main`) are defined in one task and reused consistently downstream. Artifact bases (`${study_tag}_15_scaffold`, `${study_tag}_15_<chr>_<chunk>`) and the override hooks (`GPG`, `LD_SHIP_CMD`, `LD_PREPARE_CMD`) match across the wrapper and every test. Test helper names (`_call`, `_gpg_wrapper`, `_ship_ok`, `_make_chunk`, `_prepare_stub`, `_process_chr_env`) are defined before reuse. ✓
