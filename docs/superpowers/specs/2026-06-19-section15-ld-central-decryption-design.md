# Section 15 (LD) — central decryption & reassembly design

**Date:** 2026-06-19
**Status:** Approved (design); implementation pending.
**Companion spec:** `2026-06-19-section15-ld-gpg-upload-design.md` (cohort-side
encryption + upload). This spec implements the "Central decryption (deferred)"
section of that document.

## Goal

Give central a tool to turn a cohort's uploaded, GPG-encrypted section-15 (LD)
archives back into the `cohort_stats/` directory tree that
`ld_aggregate.py accumulate` consumes — verifying integrity on the way — so a
decrypted cohort can be folded into the pooled LD precursor.

This spec covers the **central-side decryption + reassembly only**. It assumes the
`.aes`/`.md5sum` archives are already on local disk (pulled via existing sftp
tooling); decryption performs **no network I/O**, mirroring `ld_encrypt_cohort.sh`.

## Background / constraints

- **Cohort-side encryption** (`resources/genetics/ld_encrypt_cohort.sh`) produces,
  per cohort, in a flat upload dir:
  - `<study>_15_scaffold.tgz.aes` + `<study>_15_scaffold.md5sum` — a tar of the
    small scaffold files, currently `manifest.json variants.tsv.gz D.npy B.npy`.
  - one `<study>_15_chr<C>_chunk_<N>.tgz.aes` + `.md5sum` per `A_blocks` chunk —
    each tar built with `-C <cohort_stats>/A_blocks` so its contents are
    `chr<C>/chunk_<N>/…`.
  - Each `.md5sum` is `md5sum` of the **plaintext `.tgz`** (transport integrity).
- **Encryption method:** symmetric `gpg --symmetric --cipher-algo AES256`, single
  consortium passphrase shared with developers out-of-band. Decryption uses the
  same passphrase via gpg-agent cache / single prompt — nothing per-cohort.
- **Two independent integrity layers** (both must survive the round-trip):
  - `.md5sum` — md5 of the plaintext tar; checked at decrypt time.
  - `checksums.json` — blake2b over the plaintext artefacts
    (`resources/genetics/ld_checksums.py`); re-verified by `accumulate` before a
    cohort is folded in and its raw data deleted.
- **The checksums gap:** the scaffold tar does **not** currently include
  `checksums.json`, so a tree rebuilt purely from the `.aes` archives would have
  nothing for `accumulate`'s blake2b verification to check. This spec closes that
  gap (see "Cohort-side change").
- **Central accumulate interface** (`15b-ld_aggregate_stats.sh` →
  `ld_aggregate_stats.py` → `ld_aggregate.accumulate`):
  - `--cohort-dir ${ld_prepare_dir}` is the **single** directory input; it must
    contain `manifest.json` at its top, plus `variants.tsv.gz`, `B.npy`, `D.npy`,
    `A_blocks/chr*/chunk_*/`, and `checksums.json`.
  - `study_name` is **derived from `manifest.json["study_name"]`**
    (`ld_aggregate.py:440`); it is *not* a CLI argument.
  - `${ld_prepare_dir}` is a single config var, and central accumulates cohorts
    **one at a time**.

## Interface alignment (verified against `accumulate`)

| `ld_decrypt_cohort.sh` arg | Counterpart in `accumulate`                                  |
|----------------------------|--------------------------------------------------------------|
| `output_dir`               | **Exact match** for `--cohort-dir` / `${ld_prepare_dir}` — the rebuilt `cohort_stats/` tree. |
| `study_name`               | **Not an accumulate arg.** Decrypt uses it only to select a cohort's `<study>_15_*` archives in `input_dir`; accumulate reads `study_name` from the decrypted manifest. A new cross-check (below) asserts the two agree. |
| `input_dir`                | **No counterpart.** Decrypt is the pure local transform that sits *before* accumulate. |

Operator loop: *decrypt cohort A → `accumulate --cohort-dir <A>` → decrypt cohort
B → `accumulate --cohort-dir <B>` → … → finalise.* Each decrypt's `output_dir` is
the dir that run's `--cohort-dir` points at.

## Design

### New script: `resources/genetics/ld_decrypt_cohort.sh`

Exact structural mirror of `ld_encrypt_cohort.sh` (same `set -euo pipefail`, same
`GPG="${GPG:-gpg}"` test override, same logging prefix style).

```
Usage: ld_decrypt_cohort.sh <input_dir> <output_dir> <study_name>
  input_dir   dir holding <study>_15_*.tgz.aes + .md5sum (left untouched)
  output_dir  where the cohort_stats/ tree is rebuilt (== accumulate --cohort-dir)
  study_name  selects this cohort's archives by filename prefix
```

### Per-archive restore (mirror of `stage_archive`)

A `restore_archive <base>` helper, inverting `stage_archive`:

1. `"${GPG}" --output "${work}/${base}.tgz" --decrypt "${input_dir}/${base}.tgz.aes"`
   — decrypt into a staging path under `output_dir/.staging/` (the input `.aes`
   is never modified or removed).
2. **Verify md5 before untar:** `( cd "${work}" && md5sum -c "${input_dir}/${base}.md5sum" )`.
   The `.md5sum` references the bare filename `<base>.tgz`, so `md5sum -c` must run
   with cwd = the dir containing the decrypted `.tgz`. Hard-fail on mismatch.
3. **Untar atomically:** extract into a temp `*.partial` dir, then `mv` into the
   final location — so a present target dir reliably means a *complete* extraction
   (robust resume after an interruption mid-untar).
4. `rm -f "${work}/${base}.tgz"` — bounds peak disk to ~one chunk's tar.

Two extraction targets, matching how the tars were built:
- **Scaffold** (`<study>_15_scaffold`): untar into `output_dir` → drops
  `manifest.json`, `variants.tsv.gz`, `D.npy`, `B.npy`, `checksums.json` at top.
- **Chunks** (`<study>_15_chr<C>_chunk_<N>`): untar into `output_dir/A_blocks/`
  (created with `mkdir -p`) → recreates `chr<C>/chunk_<N>/…`.

### Order, enumeration, resume, guards

1. **Scaffold first.** Restore the scaffold, then read `output_dir/manifest.json`.
2. **Study-name cross-check (new):** assert
   `manifest.json["study_name"] == <study_name>`; hard-fail on mismatch (catches a
   mislabeled or cross-contaminated upload before it is folded in under the wrong
   identity). Use `jq` if available; otherwise a minimal grep/python fallback.
3. **Enumerate chunks** by globbing
   `input_dir/<study>_15_chr*_chunk_*.tgz.aes` (nullglob).
4. **Resume:** skip the scaffold if `output_dir/manifest.json` already exists; skip
   a chunk if its `output_dir/A_blocks/chr<C>/chunk_<N>/` dir already exists.
5. **Guards (mirror encrypt):** fail loudly if zero chunk archives are found; warn
   (best-effort `jq`) if the found chunk count ≠ the manifest's declared
   `[.A_blocks.chromosomes[].n_chunks] | add`.
6. Log a final `[ld_decrypt] done: <n> chunk archive(s) + scaffold in <output_dir>`.

### Cohort-side change (closes the checksums gap)

In `ld_encrypt_cohort.sh`, add `checksums.json` to the scaffold tar (line 47):

```
tar czf "${output_dir}/${scaffold}.tgz" -C "${cohort_stats_dir}" \
    manifest.json variants.tsv.gz D.npy B.npy checksums.json
```

15a writes `checksums.json` as its final step (`ld_checksums.write_cohort_checksums`),
so it is present at encrypt time. This is the only cohort-side edit; it keeps the
blake2b layer intact end-to-end so `accumulate` can verify a decrypted cohort with
`verify_checksums=True` (the default) unchanged.

### Integrity, end to end

Both layers preserved and independent:
- `.md5sum` (md5 of plaintext tar) — checked by `ld_decrypt_cohort.sh` at decrypt,
  before any untar.
- `checksums.json` (blake2b over artefacts) — checked by `accumulate` as it does
  today, now against a `checksums.json` that actually exists in the rebuilt tree.

### Error handling

- `set -euo pipefail`; any `gpg`, `md5sum -c`, or `tar` failure aborts the run.
- `mkdir -p` for `output_dir`, `output_dir/A_blocks`, and the staging dir is
  idempotent.
- Resume guards make re-running safe after an interrupted decrypt.
- Staging `.tgz` files are removed on success; the `.staging/` dir is cleaned up at
  the end (and `*.partial` dirs from a crashed untar are removed/ignored on
  re-run before re-extracting).

## Testing

New `tests/genetics/test_ld_decrypt_cohort.py`, plus an update to
`test_ld_encrypt_cohort.py`:

- **Round-trip (fake `gpg`):** a fake symmetric `gpg` (copy `.tgz`→`.aes` on
  encrypt, copy `.aes`→`.tgz` on decrypt, both via the `GPG` env override) →
  encrypt a fixture `cohort_stats/`, decrypt it, assert the rebuilt tree
  byte-matches the original (scaffold files + every `A_blocks` chunk),
  `checksums.json` is present, no leftover `.tgz`, and the input `.aes` files are
  untouched.
- **Resume:** pre-create one chunk's target dir, re-run, assert that chunk is
  skipped and others still restored.
- **md5 integrity fail:** corrupt a `.md5sum` (or the encrypted payload), assert a
  hard failure *before* any untar of that archive.
- **Study-name cross-check:** decrypt with a `study_name` that disagrees with the
  scaffold manifest's `study_name`, assert a hard failure.
- **End-to-end with accumulate:** decrypt → `ld_aggregate.accumulate(... ,
  verify_checksums=True)` succeeds (blake2b layer present).
- **Update `test_ld_encrypt_cohort.py`:** assert `checksums.json` is included in
  the scaffold archive.

## Out of scope (YAGNI)

- Downloading `.aes` from SFTP (assumed already local; operator uses existing sftp
  tooling).
- Public-key GPG (symmetric chosen for consortium consistency).
- Parallel/multi-cohort orchestration (central accumulates one cohort at a time;
  the operator loop drives the sequence).
- Regenerating `checksums.json` on central (would be circular; we transport the
  cohort's own `checksums.json` instead).
- Deleting `.aes` inputs after decrypt (kept so a later failed `accumulate` can be
  retried without re-download).
