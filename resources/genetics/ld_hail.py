"""Hail-backed genotype IO for section 15.

This module owns every code path that touches PLINK ``.bed`` entries. The
pure-Python helpers in :mod:`ld_qc` (variant index, sample alignment,
covariate matrix) run before this module is invoked so that fast unit
tests stay in pure Python.

Sex factor source-of-truth note:
    ``hl.import_plink`` populates ``mt.is_female`` from PLINK FAM
    column 5. Section 15 does NOT consume that field. ``Sex_factor`` is
    read from the covariates file via :func:`build_covariate_matrix` and
    is the only sex coding the panel uses. A future reader must not
    double-recode by also pulling ``mt.is_female``.

Allele convention:
    PLINK BIM column 5 is A1, column 6 is A2. We map ``alt = A1`` and
    ``ref = A2``. Hail's ``import_plink`` exposes ``alleles = [A2, A1]``
    so ``alleles[0] = ref`` and ``alleles[1] = alt`` which matches the
    canonical ``chr:pos:ref:alt`` variant_id used in :mod:`ld_qc`.
"""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, TypedDict

import hail as hl
import numpy as np
from hail.linalg import BlockMatrix

if TYPE_CHECKING:
    from ld_qc import VariantIndex


class VariantDiagnostics(TypedDict):
    variant_id: str
    n_nonmissing: int
    n_imputed: int
    genotype_mean: float


DEFAULT_N_PARTITIONS = 32
HAIL_VERSION = hl.version()

_INITIALIZED = False


def init_hail(
    log_file: str | Path,
    n_partitions: int = DEFAULT_N_PARTITIONS,
) -> dict[str, str]:
    """Initialise Hail once per process and return runtime metadata.

    Subsequent calls are no-ops and return the cached metadata so the
    cohort entry-point can call this freely without tracking state.
    """
    global _INITIALIZED

    if not _INITIALIZED:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        hl.init(log=str(log_path), quiet=True)
        hl.default_reference("GRCh37")
        _INITIALIZED = True

    return {
        "hail_version": HAIL_VERSION,
        "n_partitions_default": str(n_partitions),
    }


def load_genotype_matrixtable(
    bfile: str | Path,
    chromosome: str | None,
    final_samples: list[str],
    variant_index: "VariantIndex",
    n_partitions: int = DEFAULT_N_PARTITIONS,
) -> hl.MatrixTable:
    """Load a PLINK bfile as a Hail MatrixTable in canonical row/col order.

    The MatrixTable rows are filtered to the canonical ``variant_id`` set
    in ``variant_index`` and ordered to match it. The columns are
    filtered and reordered to ``final_samples`` exactly. Hard-fails if
    any expected variant_id or sample is missing from the bfile.
    """
    if not _INITIALIZED:
        raise RuntimeError(
            "init_hail must be called before load_genotype_matrixtable"
        )

    bfile_str = str(bfile)
    mt = hl.import_plink(
        bed=f"{bfile_str}.bed",
        bim=f"{bfile_str}.bim",
        fam=f"{bfile_str}.fam",
        reference_genome="GRCh37",
        skip_invalid_loci=False,
        n_partitions=n_partitions,
    )

    if chromosome is not None:
        mt = mt.filter_rows(mt.locus.contig == chromosome)

    expected_records = variant_index["variants"]
    expected_order = [v["variant_id"] for v in expected_records]
    expected_set = set(expected_order)

    mt = mt.annotate_rows(
        variant_id=hl.delimit(
            [
                mt.locus.contig,
                hl.str(mt.locus.position),
                mt.alleles[0],
                mt.alleles[1],
            ],
            ":",
        )
    )
    mt = mt.filter_rows(hl.literal(expected_set).contains(mt.variant_id))

    fam_iids = [str(s) for s in mt.s.collect()]
    iid_to_index = {iid: i for i, iid in enumerate(fam_iids)}
    missing_samples = [s for s in final_samples if s not in iid_to_index]
    if missing_samples:
        raise ValueError(
            f"final_samples missing from PLINK column set "
            f"({len(missing_samples)} samples; first: {missing_samples[:5]})"
        )
    column_order = [iid_to_index[s] for s in final_samples]
    mt = mt.choose_cols(column_order)

    actual_order = mt.variant_id.collect()
    if actual_order != expected_order:
        if len(actual_order) != len(expected_order):
            detail = (
                f"row count mismatch: Hail returned {len(actual_order)}, "
                f"variant_index expects {len(expected_order)}"
            )
        else:
            first_diff = next(
                (
                    i
                    for i, (a, b) in enumerate(zip(actual_order, expected_order))
                    if a != b
                ),
                None,
            )
            detail = f"first mismatch at row index {first_diff}"
        raise ValueError(
            f"Hail post-filter locus order does not match variant_index "
            f"({detail})"
        )

    return mt


def prepare_for_cross_products(
    mt: hl.MatrixTable,
) -> tuple[hl.MatrixTable, list[VariantDiagnostics]]:
    """Annotate per-variant diagnostics, mean-impute genotype dosages.

    Returns the MatrixTable with a new ``GT_dosage`` entry field (float64;
    missing calls replaced with the per-variant mean alt-allele count) and
    a list of ``VariantDiagnostics`` in canonical row order. Hard-fails if
    any variant has zero non-missing calls, since mean imputation would
    produce NaN dosages and corrupt downstream ``B``/``A`` cross-products.
    """
    n_cols = mt.count_cols()
    mt = mt.annotate_rows(
        n_nonmissing=hl.agg.count_where(hl.is_defined(mt.GT)),
        genotype_mean=hl.agg.mean(mt.GT.n_alt_alleles()),
    )
    mt = mt.annotate_rows(n_imputed=n_cols - mt.n_nonmissing)

    rows = (
        mt.rows()
        .select("variant_id", "n_nonmissing", "n_imputed", "genotype_mean")
        .collect()
    )
    all_missing = [r.variant_id for r in rows if r.n_nonmissing == 0]
    if all_missing:
        raise ValueError(
            f"{len(all_missing)} variants are entirely missing across all "
            f"{n_cols} samples; first 5: {all_missing[:5]}"
        )

    diagnostics: list[VariantDiagnostics] = [
        {
            "variant_id": r.variant_id,
            "n_nonmissing": int(r.n_nonmissing),
            "n_imputed": int(r.n_imputed),
            "genotype_mean": float(r.genotype_mean),
        }
        for r in rows
    ]

    mt = mt.annotate_entries(
        GT_dosage=hl.coalesce(
            hl.float64(mt.GT.n_alt_alleles()),
            mt.genotype_mean,
        )
    )
    return mt, diagnostics


def compute_b_block(
    mt: hl.MatrixTable,
    covariate_matrix: np.ndarray,
) -> np.ndarray:
    """Compute ``B_k = X^T C`` (variants × covariates) as a NumPy array.

    The MatrixTable must already carry a ``GT_dosage`` entry field (call
    :func:`prepare_for_cross_products` first). The rows of
    ``covariate_matrix`` must align one-to-one with the columns of ``mt``
    in the same order (i.e. ``covariate_matrix[i]`` corresponds to the
    sample in ``mt`` column ``i``).

    Returns a contiguous ``(n_variants, n_covariates)`` float64 array.
    The intermediate Hail ``BlockMatrix`` for the genotype matrix stays
    partitioned and is never materialised on the driver.
    """
    n_samples = mt.count_cols()
    if covariate_matrix.shape[0] != n_samples:
        raise ValueError(
            f"covariate_matrix has {covariate_matrix.shape[0]} rows but "
            f"MatrixTable has {n_samples} samples; they must match exactly"
        )
    if covariate_matrix.dtype != np.float64:
        raise ValueError(
            f"covariate_matrix must be float64; got {covariate_matrix.dtype}"
        )

    x_h = BlockMatrix.from_entry_expr(mt.GT_dosage)
    c_bm = BlockMatrix.from_numpy(np.ascontiguousarray(covariate_matrix))
    return (x_h @ c_bm).to_numpy()


DEFAULT_LD_RADIUS_BP = 10_000_000
DEFAULT_A_BLOCK_SIZE = 4096


def compute_a_block_banded(
    mt: hl.MatrixTable,
    variant_index: "VariantIndex",
    out_dir: str | Path,
    radius_bp: int = DEFAULT_LD_RADIUS_BP,
    block_size: int = DEFAULT_A_BLOCK_SIZE,
) -> dict:
    """Compute and write the upper-triangular banded ``A_k = X X^T`` per chromosome.

    For each chromosome present in ``variant_index``, the MatrixTable is
    filtered to that chromosome, ``X`` is built as a Hail BlockMatrix from
    the imputed ``GT_dosage`` entries, the cross-product ``X @ X.T`` is
    sparsified to the upper triangle within ``radius_bp`` physical
    distance, and the result is written to ``{out_dir}/chr<C>/`` as a
    Hail BlockMatrix directory. The dense ``p × p`` matrix is never
    materialised.

    Returns a manifest dict describing the written artefacts.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    by_chrom: dict[str, list[int]] = {}
    chrom_order: list[str] = []
    for variant in variant_index["variants"]:
        chrom = variant["chr"]
        if chrom not in by_chrom:
            by_chrom[chrom] = []
            chrom_order.append(chrom)
        by_chrom[chrom].append(variant["pos"])

    chromosomes_meta: dict[str, dict] = {}
    for chrom in chrom_order:
        positions = np.asarray(by_chrom[chrom], dtype=np.int64)
        n_variants = positions.size

        j_max_inclusive = (
            np.searchsorted(positions, positions + radius_bp, side="right") - 1
        )
        idx = np.arange(n_variants, dtype=np.int64)
        max_idx_distance = int((j_max_inclusive - idx).max())

        mt_chr = mt.filter_rows(mt.locus.contig == chrom)
        x_h = BlockMatrix.from_entry_expr(mt_chr.GT_dosage, block_size=block_size)
        a = x_h @ x_h.T
        a_band = a.sparsify_band(
            lower=0, upper=max_idx_distance, blocks_only=True
        )
        chr_dir = out_path / f"chr{chrom}"
        a_band.write(str(chr_dir), overwrite=True)

        chromosomes_meta[chrom] = {
            "directory": f"{out_path.name}/chr{chrom}",
            "n_variants": int(n_variants),
            "max_idx_distance_in_band": max_idx_distance,
            "block_size": block_size,
        }

    return {
        "directory": out_path.name,
        "format": "hail-blockmatrix",
        "radius_bp": radius_bp,
        "block_size": block_size,
        "chromosomes": chromosomes_meta,
    }
