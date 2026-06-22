# Section 15 (LD) — cohort run-upload-delete wrapper design

**Date:** 2026-06-22
**Status:** Approved (cohort-side wrapper + Imperial transport + wiki). Central
reassembly/decryption is an explicit follow-up, NOT in this scope.

## Goal

Give a GoDMC cohort a single, resumable command to produce, ship, and locally
reclaim their section-15 (LD) cohort outputs, without ever needing the ~2 TB of
free disk a whole-genome run would otherwise require. The deliverable is both the
wrapper script and a cohort-facing wiki runbook describing how to run it.

## Why section 15 is different from every other module

Every other pipeline module produces small outputs, so the cohort convention is
"produce everything → `check_upload.sh <n> check` → `check_upload.sh <n> upload`".
Section 15 cannot follow this. Its genome-wide `A_blocks` output is ~1.3 TB per
cohort (chr1 alone ≈ 103 GB; genome-wide extrapolation ≈ 1.3 TB / ~33 h — see
`LD_IMPLEMENTATION_PROGRESS.md`). **The disk pressure is during `15a` production,
not at upload.** By the time a whole-genome `15a` has finished, the 1.3 TB is
already on disk. Therefore production and upload-delete must be *interleaved per
chromosome*; they cannot be two separate whole-genome phases.

## What is big vs. small

- **Scaffold** (`D.npy`, `B.npy`, `variants.tsv.gz`, `manifest.json`,
  `checksums.json`, `qc_report.txt`): small. Intercept-only `B` is
  ~30M variants × 1 float64 ≈ 240 MB genome-wide; `D` is trivial; `variants.tsv.gz`
  is small. The whole-genome scaffold is well under 1 GB.
- **`A_blocks`**: effectively all of the ~1.3 TB.

Consequence: only `A_blocks` ever needs per-chromosome deletion. Per-chromosome
scaffolds are cheap and can accumulate locally as a record of what was shipped.

## Approach (chosen)

A new cohort-side orchestrator script, **`15c-ld_run_upload.sh`**, is the single
command a cohort operator runs. It does NOT replace `15a-ld_prepare_stats.sh` — it
calls the existing per-chromosome compute worker once per chromosome, pointing each
chromosome at its own output directory, and interleaves encrypt → upload → verify →
delete so peak disk stays at roughly one chromosome's `A_blocks` (~100 GB) rather
than the full ~1.3 TB.

Rejected alternative — a **documented manual per-chromosome bash loop** using the
existing `15a` + `check_upload.sh 15 upload` + `rm`: it puts `rm -rf` of
unrecoverable data in the operator's hands, is not resumable without bookkeeping,
and still requires the same underlying fixes (scaffold overwrite, encrypt-helper
scaffold name collision). Same fix cost, worse safety.

Rejected alternative — **folding the loop into `check_upload.sh`'s section-15
block**: `check_upload.sh` operates on already-produced outputs, so it cannot bound
production-time disk; and section 15 ships to Imperial, not the Bristol SFTP that
`check_upload.sh` targets.

## Design

### Components and responsibilities

- **`15a-ld_prepare_stats.py`** (unchanged) — per-chromosome compute worker. Already
  accepts `--chromosome <C>` and `--output-dir <dir>`. The wrapper drives it once
  per chromosome with `--output-dir results/15/cohort_stats/chr<C>/`. Pointing each
  chromosome at its own directory dissolves the scaffold-overwrite problem with no
  code change to the worker.
- **`resources/genetics/ld_encrypt_cohort.sh`** (small tweak) — already stages +
  symmetric-GPG-encrypts a `cohort_stats` directory into per-chunk `.tgz.aes` +
  `.md5sum`, with a resume guard. The single fixed scaffold name
  (`${study_name}_15_scaffold`) collides across chromosomes, and the resume guard
  would skip chromosomes 2–22's scaffold. Fix: include the chromosome in the
  artifact names. Cleanest: the wrapper invokes the helper with a per-chromosome
  study tag (`${study_name}_chr<C>`) so every artifact — scaffold and chunks — is
  chromosome-namespaced, with no change to the helper's internal logic.
- **`15c-ld_run_upload.sh`** (new) — the cohort entry point. Sources config /
  `resources/parameters`, then runs the per-chromosome loop below.
- **Wiki** — `Run-federated-LD-reference-panel.md` rewritten as a cohort operator
  runbook (see "Wiki" section).

### The per-chromosome loop

For each chromosome `C` in `${ld_chromosomes}` (default `1..22`, overridable for
pilots):

1. **Skip if done.** If `results/15/cohort_stats/chr<C>/.uploaded` exists, skip the
   chromosome. This makes the entire loop resumable after any interruption.
2. **Compute.** Run `15a` for chromosome `C` into `cohort_stats/chr<C>/`, producing
   that chromosome's scaffold + `A_blocks/chr<C>/chunk_*/` + `checksums.json`.
3. **Check.** Run the section-15 checks scoped to this chromosome (outputs exist,
   QC sane). Abort the chromosome on failure without deleting anything.
4. **Encrypt + upload + delete, streamed per chunk.** For each `A_blocks` chunk:
   encrypt to `.tgz.aes` (the helper `rm`s the intermediate `.tgz`), rsync the
   `.aes` + `.md5sum` to Imperial, verify the transfer, then `rm` both the chunk
   directory and its local `.aes`. Streaming per chunk (not per whole chromosome)
   keeps the encrypted copy from doubling peak disk — extra overhead stays at
   ~one chunk (~130 MB) on top of the chromosome's `A_blocks`.
5. **Ship scaffold + mark + reclaim.** Encrypt and upload the small scaffold
   archive, verify, then write `chr<C>/.uploaded` and `rm -rf chr<C>/A_blocks`.
   Keep the tiny scaffold directory locally as a record of what was shipped.

**Peak disk ≈ one chromosome's `A_blocks` (~100 GB for chr1) + ~one chunk.**

### Safety: delete only after verified upload

The single non-negotiable rule: **nothing is deleted until its upload is verified.**
`15a` already emits `checksums.json` (stdlib blake2b). Verification compares the
shipped artifact against its recorded checksum; rsync runs with `--checksum` for
byte-identical transfer. Remote-side recompute is a later hardening once the
concrete Imperial endpoint is known. An interrupted run leaves the source intact
(no sentinel written ⇒ chromosome re-runs cleanly on the next invocation).

### Transport (Imperial)

rsync-over-SSH, isolated in a single function so the concrete endpoint can be
swapped without touching the loop. Flags: `--partial --append` (resumable mid-chunk
over a flaky link) and `--checksum` (byte-identity). Endpoint configuration lives in
**new** config vars so the Bristol `sftp_*` vars and `check_upload.sh` are untouched:

- `imperial_host`, `imperial_user`, `imperial_path`, `imperial_key`.

The transport is intentionally pluggable: the Imperial endpoint is not finalised
(probably rsync/SSH). The single-function isolation is the deliberate hedge.

### Encryption

Unchanged from the existing convention and the section-15 GPG upload spec:
symmetric `gpg --symmetric --cipher-algo AES256`, passphrase shared with developers
out-of-band, gpg-agent caches it within the session (one prompt in practice).
Per-chunk archives. The only change is chromosome-namespaced artifact names.

### Config / parameters

New section-15 cohort vars (defaults in `resources/parameters`):

- `ld_chromosomes` — chromosome list, default `1..22`. Override for pilots
  (e.g. `22`).
- `imperial_host`, `imperial_user`, `imperial_path`, `imperial_key` — Imperial
  rsync endpoint.

Existing section-15 vars (`ld_chromosome`, `ld_hail_*`, `ld_a_*`,
`ld_prepare_dir`, etc.) are reused; the wrapper sets per-chromosome `--output-dir`
itself rather than relying on a single `ld_prepare_dir`.

## Wiki

Rewrite `godmc_phase2.wiki/Run-federated-LD-reference-panel.md` from developer notes
into a cohort operator runbook:

- **Prerequisites:** sections `00`, `01`, `02a`, `03a`; the Hail env
  (`resources/genetics/hail_env.yml`); **~150 GB free disk (not 2 TB)**.
- **The one command:** `./15c-ld_run_upload.sh <config>`.
- **What it does:** plain-English description of the per-chromosome
  run → check → encrypt → upload → verify → delete loop, and why section 15 works
  this way (disk).
- **Configuration:** the new Imperial vars, `ld_chromosomes`, and passphrase
  handling.
- **Resumability:** safe to re-run after interruption; already-shipped chromosomes
  are skipped.
- **What is deleted and when:** `A_blocks` per chromosome only after its upload is
  verified; scaffolds kept locally as a record.
- **Expectations:** chr1 ≈ 2.5 h / ~100 GB peak; genome-wide ≈ 33 h streamed.
- **Progress / recovery:** how to see which chromosomes are done (sentinels), how to
  resume.
- Update the stale `Scripts status` ("In development") and `Data upload method`
  ("not yet implemented") headers.

## Testing

- **Wrapper parse:** `bash -n 15c-ld_run_upload.sh`.
- **Loop logic (mocked compute + transport):** with `15a` and the rsync step
  stubbed, assert the loop (a) skips chromosomes with a `.uploaded` sentinel,
  (b) does NOT delete `A_blocks` when the upload/verify step fails, (c) writes the
  sentinel and deletes `A_blocks` only after a successful verified upload,
  (d) namespaces artifacts per chromosome (no scaffold collision across chrs).
- **Encrypt-helper per-chr naming:** running the helper with a `${study}_chr<C>`
  tag yields chromosome-namespaced scaffold + chunk artifacts; two chromosomes do
  not collide and neither is skipped by the resume guard.
- **Round-trip sanity (manual, pilot):** on the chr22/chr1 pilot output, run the
  wrapper against a local rsync target (`rsync` to a local path), confirm the
  decrypt of one chunk reproduces the source byte-for-byte, and confirm
  `A_blocks` is gone locally afterwards while the sentinel + scaffold remain.

## Out of scope (follow-up)

- **Central reassembly / decryption.** Per-chromosome output means the central side
  receives per-chr scaffolds; `ld_decrypt_cohort.sh` / 15b must concatenate
  `variants`/`B` rows across chromosomes and validate-then-take a single `D`
  (covariate-only, identical across chromosomes). This is real central-side work
  and lands as a **separate spec/plan**. It must be complete before the first real
  cohort runs end-to-end, but it does not block the wrapper or the wiki.
- Remote-side checksum recompute (hardening once the Imperial endpoint is fixed).
- Non-rsync transports (the transport function is the swap point if Imperial
  finalises on something else).
- Parallel/multi-chromosome concurrency (the loop is deliberately serial to bound
  disk).
