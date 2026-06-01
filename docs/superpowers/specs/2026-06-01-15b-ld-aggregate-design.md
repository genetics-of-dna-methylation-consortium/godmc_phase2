# Section 15b — Central LD Aggregation Design

Date: 2026-06-01
Status: Approved design, pre-implementation
Companion documents: `LD_IMPLEMENTATION_PLAN.md`, `LD_IMPLEMENTATION_PROGRESS.md`

## Aim

Replace the pre-Hail 15b scaffold (`resources/genetics/ld_aggregate_stats.py`) with a
production central aggregator that turns per-cohort sufficient statistics (`A_k`, `B_k`,
`D_k`, exported by the feature-complete 15a) into a pooled covariate-adjusted LD
correlation panel `R`, without ever centrally storing individual-level genotypes.

The target estimand is unchanged from the plan: **full pooled covariate-adjusted LD from
the stacked cohort genotype matrix**, restricted to the variant set common to all cohorts.

## Context and inputs

15a is feature-complete and, per cohort, writes under its output directory:

- `manifest.json` — covariate schema, `genome_build`, `radius_bp`, `A_blocks.block_size`,
  variants schema version, counts.
- `variants.tsv.gz` — canonical `chr:pos:ref:alt` order with per-variant `n_nonmissing`,
  `n_imputed`, `genotype_mean`.
- `D.npy` — dense `D_k = C_k^T C_k`, 3×3, float64.
- `B.npy` — dense `B_k = X_k^T C_k`, (variants × 3), float64, in `variants.tsv.gz` order.
- `A_blocks/chr<C>/chunk_<N>/` — Hail BlockMatrix dirs, upper-triangular row-interval
  `A_k = X_k X_k^T` within a 1 Mb radius, float64.

The covariate matrix columns are `[intercept, Age_numeric, Sex_factor]` (`Sex_factor`
recoded `M=1, F=2`). `A`, `B`, `D` are raw-genotype-scale (uncentred, unstandardised);
all centring/scaling happens centrally.

## Decisions (locked during brainstorming)

1. **Pooled intersection on the canonical key.** Because every cohort imputes against the
   same panel, `chr:pos:ref:alt` keys are directly comparable with no allele-flip or
   liftover rescue. But per-cohort section-02 QC (INFO, MAF, missingness, HWE) drops
   different variants, so variant sets overlap heavily but are not byte-identical.
   Therefore the pooled set is the **intersection** of variants present in *all* cohorts,
   and each cohort is re-indexed onto a frozen pooled index. This preserves the full-pooled
   estimand: every retained entry sums the identical set of cohorts.

2. **Chunk-streaming NumPy/SciPy engine.** All aggregation, covariate adjustment, and `R`
   conversion run as bounded per-chromosome row-chunks (mirroring 15a's writer), not as
   native Hail BlockMatrix arithmetic. The central math (`A_adj = A − B D⁻¹ Bᵀ`, then `R`)
   is fundamentally entrywise within the stored 1 Mb window; the `B D⁻¹ Bᵀ` correction is
   rank-3, so forming it in BlockMatrix space risks a dense p×p materialisation, while a
   bounded NumPy block computes it exactly within the window. Hail is still used to *read*
   cohort `A_blocks` back (windowed) and to *write* the final `R_blocks`.

3. **`R` output format.** Float32 Hail BlockMatrix chunks, identical layout to 15a's
   `A_blocks` (`R_blocks/chr<C>/chunk_<N>/`, upper-triangular row-interval). The plan
   permits downcasting the final `R` to float32; everything upstream of `R` stays float64.

4. **Defer all filtering to finalisation.** Pooled MAF/MAC and unstable-diagonal filters
   are applied centrally after aggregation (never per-cohort). No global ridge in release 1;
   PSD is enforced only on the small submatrices extracted later for fine-mapping. This
   keeps the running precursor as pure additive sufficient statistics.

5. **Two-phase incremental design.** 15b is `accumulate` (add one cohort to a mutable
   precursor) + `finalise` (emit an immutable versioned panel). This bounds central peak
   storage to roughly one precursor + one incoming cohort regardless of cohort count `K`,
   and lets new cohorts join without reprocessing prior ones.

6. **No per-file checksums.** Nothing else in the pipeline hashes files; 15b validates
   cohorts via `manifest.json` contract fields and identifies cohorts by `study_name`. The
   previously deferred "add sha256 checksums" item is resolved as won't-do.

## Architecture

15b has two entry points, both wrapped by `15b-ld_aggregate_stats.sh` (dispatched by a
`--mode {accumulate,finalise}` flag):

- **`accumulate`** — adds exactly one cohort's `A`/`B`/`D` into the running precursor.
  Re-runnable; refuses an already-present `study_name` unless `--force`.
- **`finalise`** — reads the precursor, resolves the intersection, filters, computes
  `A_adj` and `R`, and writes an immutable versioned panel. Re-runnable; each run emits a
  new panel version.

### Precursor state (mutable, on disk)

Under `${ld_precursor_dir}` (default `${ld_aggregate_dir}/precursor/`):

| File | Content |
|---|---|
| `precursor_manifest.json` | panel-spec version, covariate-schema id, `radius_bp`, `block_size` contract, `genome_build`, accumulated-cohort list (`study_name`, timestamp, variant count), total cohort count `K` |
| `D.npy` | running Σ `D_k`, 3×3 float64 |
| `variants.parquet` | append-only master union table: `stable_id` (monotonic int, never reused), `chr/pos/ref/alt`, `membership_count`, the three accumulated `B` columns (`b_intercept`, `b_age`, `b_sex`), accumulated `A` diagonal `a_diag`, accumulated `n_nonmissing`, `n_imputed` |
| `A_pairs/chr<C>.parquet` | off-diagonal windowed `A` sums only, COO `(pos_i, sid_i, pos_j, sid_j, value)` **sorted by `(pos_i, sid_i, pos_j, sid_j)`** (i.e. genomic-position upper triangle, `pos_i ≤ pos_j`), mergeable by linear sorted-merge-add |

Design rationale: everything *per-variant* that is additive (the `B` row, the `A`
diagonal, frequency inputs, membership count) lives in `variants.parquet`; the only large
structure is the **off-diagonal pair store**. The variant table carries a `stable_id`
assigned at first sighting (append-only) so the union can grow **without re-indexing
existing per-variant data**. The pair store is *sorted by genomic position*
(`pos_i, sid_i, pos_j, sid_j`); genomic position is stable as the union grows, equals the
pooled ordering, and is the order in which a cohort's `A_blocks` are read back chunk by
chunk — so `accumulate` is a linear sorted-merge and `finalise` streams the pair store in
pooled-row order without ever sorting a ~10¹⁰-row file. The `sid` columns disambiguate
multiple variants at one position and map pairs back to the variant table; mapping to
contiguous pooled indices happens lazily inside `finalise`.

### Storage lifecycle

Peak central footprint ≈ precursor (~one cohort's windowed `A`) + one incoming cohort,
independent of `K`. After `accumulate` records the cohort in `precursor_manifest.cohorts`,
the cohort's raw `A_blocks` upload can be deleted.

## `accumulate` phase

Inputs: one 15a cohort directory + the precursor directory (created on first call).

1. **Validate contract** (hard-fail on mismatch): covariate `schema_id` / `matrix_columns`
   / sex recode, `genome_build`, `radius_bp`, `A_blocks.block_size`, variants schema
   version. On a fresh precursor the first cohort *defines* these values; every later
   cohort must match.
2. **Idempotency guard:** if this `study_name` is already in
   `precursor_manifest.cohorts`, refuse. `--force` allows a deliberate re-run only after a
   `remove`/rebuild.
3. **Merge per-variant data** (streamed in canonical-key order): look up each cohort
   variant_id in the master table, append a new `stable_id` if unseen; add the three `B`
   columns, bump `membership_count`, add `n_nonmissing`/`n_imputed`.
4. **Merge `A`** by reading cohort `A_blocks` chunks back as NumPy in chunk order (each
   chunk is already bounded by 15a's `max_dense_gb`). Chunk row/col indices are
   chromosome-relative position order, so entries arrive sorted by `(pos_i, pos_j)`.
   Diagonal entries (`j == i`) add into `variants.parquet` `a_diag`; off-diagonal entries
   become `(pos_i, sid_i, pos_j, sid_j, value)` rows that **linear sorted-merge-add** into
   `A_pairs/chr<C>.parquet` (one streaming pass over the existing precursor file plus the
   incoming sorted stream).
5. **Add `D_k`** into precursor `D.npy`; increment `K`.
6. **Commit:** append the cohort record to `precursor_manifest.json` and write atomically
   (write-temp-then-rename). The manifest rename is the commit point; the pair-store merge
   also writes to a temp parquet and renames last, so a crash before commit leaves the
   precursor unchanged.

Concurrency/safety: `accumulate` holds a lockfile on the precursor directory; concurrent
accumulate fails fast. All writes are temp-then-rename.

## `finalise` phase

Inputs: the precursor directory. Output: `${ld_panel_dir}/panel_v<N>/` (immutable).

1. **Resolve the intersection:** keep variants with `membership_count == K` (default; a
   `--min-cohorts` param can relax this). Sort survivors by `(chr, pos, ref, alt)` into
   contiguous **pooled indices**; build the `stable_id → pooled_index` map.
2. **Assemble small pooled objects:** `B` (m×3), `D` (3×3), `A` diagonal (m-vector), all
   for the intersection.
3. **Filter (with recorded per-variant reasons):**
   - pooled ALT freq `f = B[:,0] / (2·D[0,0])`, `MAF = min(f, 1−f)`; drop
     `MAF < --maf-threshold` (default 0.01),
   - adjusted diagonal `A_adj[i,i] = a_diag[i] − B[i]·(D⁻¹ B[i]ᵀ)`; drop non-positive or
     `< --min-adj-diag` entries,
   - re-prune the pooled index after dropping; log all drops to `dropped_variants.tsv`.
4. **Pre-solve once:** check `D` rank/condition (hard-fail by default; optional
   ridge-on-`D` only if explicitly requested). Form `W = (D⁻¹ Bᵀ)ᵀ` (m×3) so the
   correction is `corr[i,j] = B[i]·W[j]`.
5. **Stream per chromosome in bounded pooled row-chunks** (same `max_dense_gb` logic as
   15a). Because `A_pairs/chr<C>.parquet` is sorted by `(pos_i, …)` and the pooled index is
   the position rank, the pair store is read sequentially in lockstep with chunk iteration
   — no global re-sort:
   - for chunk rows `[r0, r1)` and 1 Mb window cols `[r0, c1)`: take the next run of pair
     rows whose `pos_i` falls in this chunk, drop any pair with an endpoint outside the
     intersection, map endpoints to pooled indices, and scatter into a dense bounded block;
     add the diagonal,
   - subtract the rank-3 correction block `B[r0:r1] @ W[r0:c1]ᵀ` → `A_adj` block,
   - normalise `R = A_adj / sqrt(outer(A_adj_diag[r0:r1], A_adj_diag[r0:c1]))`,
   - sparsify to the upper-triangular row intervals, downcast to **float32**, write
     `R_blocks/chr<C>/chunk_<N>/` as a Hail BlockMatrix dir.
6. **Write panel artefacts:** `pooled_manifest.json` (panel version, cohorts included,
   intersection/filter counts, `D` rank/condition, thresholds, regularisation record),
   `variants.tsv.gz` (final pooled index + pooled freq + `A_adj` diagonal + diagnostics),
   `cohort_inclusion.tsv`, `dropped_variants.tsv`, `qc_report.txt`.

Key property: `A_adj` and `R` are computed entrywise within the bounded window. `D⁻¹` is
3×3 and `W` is m×3, so the correction is rank-3 and never forms a dense p×p object. The
diagonal-in-variant-table choice makes the sqrt normalisation a clean outer-product slice.

## Error handling

Hard-fails with clear messages:

- contract mismatch on `accumulate` (covariate schema / `radius_bp` / `block_size` /
  build / variants schema version),
- re-accumulating an already-present `study_name` without `--force`,
- a cohort `A_blocks` chunk whose `block_size`/`radius_bp` disagrees with the precursor,
- `finalise` with `K == 0` (empty precursor) or an empty intersection,
- singular/ill-conditioned pooled `D` (hard-fail by default; ridge-on-`D` only when
  explicitly requested),
- lockfile contention ("another accumulate is in progress").

## Validation tests

Extend the existing local-only `tests/genetics/` harness (same `hail_env`):

1. **Drop the `pytest.skip`** in
   `test_ld_aggregate_synthetic.py::test_two_cohort_a_adj_matches_stacked_residualisation`:
   run two synthetic cohorts through 15a, `accumulate` both, `finalise`, and assert pooled
   `A_adj`/`R` match the direct stacked-residualisation NumPy reference at `rtol≈1e-12`.
2. **Incremental == batch == order-independent:** accumulate A-then-B, B-then-A, and a
   single batch all yield identical precursors (modulo timestamps) and identical `R`.
3. **Intersection correctness:** a third cohort missing some variants → those variants are
   excluded; surviving pairs sum all `K` cohorts.
4. **Re-index correctness:** `stable_id → pooled_index` mapping reproduces the right
   windowed entries when the union grew out of position order.
5. **Filter behaviour:** a rare variant is dropped at the MAF floor; a non-positive
   `A_adj` diagonal is dropped with a recorded reason.
6. **Guards:** contract-mismatch and duplicate-`study_name` hard-fails; idempotency guard
   refuses a re-run.
7. **No dense p×p:** the finalise stream never allocates a block exceeding `max_dense_gb`.

## Pipeline integration

- `15b-ld_aggregate_stats.sh` dispatches `accumulate`/`finalise` via `--mode`, defaulting
  sensibly for the two-cohort pilot.
- `resources/parameters`: add `ld_precursor_dir`, `ld_panel_dir`, `ld_maf_threshold`
  (default 0.01), `ld_min_adj_diag`, `ld_min_cohorts`, and the mode variable.
- `resources/logs/check_logs.sh` / `check_results.sh`: update success markers and result
  checks for precursor and panel outputs.
- `godmc_phase2.wiki/Run-federated-LD-reference-panel.md`: document the two-phase central
  workflow.
- This realises the **two-cohort 15b pilot** named as the next step in the progress
  tracker.

## Out of scope (release 1)

- Globally projected genotype PCs in the pooled covariate schema (deferred panel spec).
- Cohort fixed effects, age×sex, age² (deferred covariates).
- Global ridge / eigenvalue-clipping of the full panel (PSD enforced only at fine-mapping
  extraction).
- `--min-cohorts < K` is supported as a parameter but the default and validated path is
  the true full intersection.
</content>
</invoke>
