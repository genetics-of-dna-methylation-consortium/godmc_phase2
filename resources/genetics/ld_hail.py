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

import gzip
import json
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

import psutil

import hail as hl
import numpy as np
from hail.linalg import BlockMatrix

if TYPE_CHECKING:
    from ld_qc import VariantIndex


DEFAULT_N_PARTITIONS = 256
HAIL_VERSION = hl.version()
AUTOSOMES = {str(c) for c in range(1, 23)}
VALID_BASES = {"A", "C", "G", "T"}
MHC_CHROMOSOME = "6"
MHC_START_BP = 28477797
MHC_END_BP = 33448354

_INITIALIZED = False


def init_hail(
    log_file: str | Path,
    n_partitions: int = DEFAULT_N_PARTITIONS,
    local_cores: int | None = None,
    driver_memory_gb: int | None = None,
    tmp_dir: str | Path | None = None,
) -> dict[str, str]:
    """Initialise Hail once per process and return runtime metadata.

    Subsequent calls are no-ops and return the cached metadata so the
    cohort entry-point can call this freely without tracking state.
    """
    global _INITIALIZED

    if not _INITIALIZED:
        log_path = Path(log_file)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        init_kwargs: dict = {"log": str(log_path), "quiet": True}
        spark_conf: dict[str, str] = {}

        if driver_memory_gb is None:
            avail_gb = psutil.virtual_memory().available / (1024**3)
            # 0.65 rather than 0.8: leaves headroom for JVM off-heap (Netty,
            # Tungsten direct buffers, metaspace, JIT) which sit outside -Xmx
            driver_memory_gb = int(avail_gb * 0.65)

        if local_cores is not None:
            if local_cores <= 0:
                raise ValueError(f"local_cores must be positive; got {local_cores}")
            init_kwargs["backend"] = "spark"
            init_kwargs["local"] = f"local[{local_cores}]"

        if driver_memory_gb is not None:
            if driver_memory_gb <= 0:
                raise ValueError(
                    f"driver_memory_gb must be positive; got {driver_memory_gb}"
                )
            spark_conf["spark.driver.memory"] = f"{driver_memory_gb}g"

        # G1GC handles long-lived large-object heap (BlockMatrix chunks) better
        # than default GC; start a collection earlier to avoid full-GC stalls.
        spark_conf["spark.driver.extraJavaOptions"] = (
            "-XX:+UseG1GC "
            "-XX:InitiatingHeapOccupancyPercent=65 "
            "-XX:MaxGCPauseMillis=500"
        )
        # Cap Spark driver history retention to prevent metadata accumulation
        # across the hundreds of stages produced by per-chromosome chunk writes.
        spark_conf["spark.ui.retainedJobs"] = "50"
        spark_conf["spark.ui.retainedStages"] = "100"
        spark_conf["spark.ui.retainedTasks"] = "1000"

        if spark_conf:
            init_kwargs["spark_conf"] = spark_conf

        if tmp_dir is not None:
            tmp_path = Path(tmp_dir)
            tmp_path.mkdir(parents=True, exist_ok=True)
            local_tmp_path = tmp_path / "local"
            local_tmp_path.mkdir(parents=True, exist_ok=True)
            init_kwargs["tmp_dir"] = str(tmp_path)
            init_kwargs["local_tmpdir"] = str(local_tmp_path)

        hl.init(**init_kwargs)
        hl.default_reference("GRCh37")
        _INITIALIZED = True

    return {
        "hail_version": HAIL_VERSION,
        "n_partitions_default": str(n_partitions),
        "local_cores": "default" if local_cores is None else str(local_cores),
        "driver_memory_gb": (
            "default" if driver_memory_gb is None else str(driver_memory_gb)
        ),
        "tmp_dir": "default" if tmp_dir is None else str(tmp_dir),
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
    autosomes = hl.literal(AUTOSOMES)
    valid_bases = hl.literal(VALID_BASES)
    row_is_valid_snp = (
        autosomes.contains(mt.locus.contig)
        & valid_bases.contains(mt.alleles[0])
        & valid_bases.contains(mt.alleles[1])
        & (mt.alleles[0] != mt.alleles[1])
    )
    row_is_outside_mhc = ~(
        (mt.locus.contig == MHC_CHROMOSOME)
        & (mt.locus.position >= MHC_START_BP)
        & (mt.locus.position <= MHC_END_BP)
    )
    mt = mt.filter_rows(row_is_valid_snp & row_is_outside_mhc)

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

    return mt


def prepare_for_cross_products(
    mt: hl.MatrixTable,
    n_samples: int,
) -> tuple[hl.MatrixTable, hl.Table]:
    """Annotate per-variant diagnostics, mean-impute genotype dosages.

    Returns the MatrixTable with a new ``GT_dosage`` entry field (float64;
    missing calls replaced with the per-variant mean alt-allele count), a
    Hail Table of per-variant diagnostics in row order. Diagnostics are
    returned as a lazy table so full-chromosome pilots do not collect or
    aggregate millions of rows on the driver before the matrix work starts.
    """
    mt = mt.annotate_rows(
        n_nonmissing=hl.agg.count_where(hl.is_defined(mt.GT)),
        genotype_mean=hl.agg.mean(mt.GT.n_alt_alleles()),
    )
    mt = mt.annotate_rows(n_imputed=n_samples - mt.n_nonmissing)

    mt = mt.annotate_entries(
        GT_dosage=hl.coalesce(
            hl.float64(mt.GT.n_alt_alleles()),
            mt.genotype_mean,
        )
    )
    rows = mt.rows()
    diagnostics = rows.select(
        chr=rows.locus.contig,
        pos=rows.locus.position,
        ref=rows.alleles[0],
        alt=rows.alleles[1],
        variant_id=rows.variant_id,
        n_nonmissing=rows.n_nonmissing,
        n_imputed=rows.n_imputed,
        genotype_mean=rows.genotype_mean,
    )
    return mt, diagnostics


def compute_b_block(
    mt: hl.MatrixTable,
    covariate_matrix: np.ndarray,
    temp_dir: str | Path,
) -> np.ndarray:
    """Compute ``B_k = X^T C`` (variants × covariates) as a NumPy array.

    The MatrixTable must already carry a ``GT_dosage`` entry field (call
    :func:`prepare_for_cross_products` first). The rows of
    ``covariate_matrix`` must align one-to-one with the columns of ``mt``
    in the same order (i.e. ``covariate_matrix[i]`` corresponds to the
    sample in ``mt`` column ``i``).

    Returns a contiguous ``(n_variants, n_covariates)`` float64 array. This
    deliberately avoids ``BlockMatrix.from_entry_expr`` because that path
    materialises a dense ``variants × samples`` genotype matrix before the
    tiny ``variants × covariates`` result is available.

    The per-variant b-vectors are written to a directory of per-partition
    TSV files via ``parallel="separate_header"`` and concatenated in
    partition order on the Python side. A single-file export would force
    Hail to collect every row to the driver, which OOMs the driver heap on
    whole-chromosome inputs.
    """
    if covariate_matrix.dtype != np.float64:
        raise ValueError(
            f"covariate_matrix must be float64; got {covariate_matrix.dtype}"
        )
    n_samples = mt.count_cols()
    if covariate_matrix.shape[0] != n_samples:
        raise ValueError(
            f"covariate_matrix has {covariate_matrix.shape[0]} rows but "
            f"MatrixTable has {n_samples} samples; they must match exactly"
        )

    temp_path = Path(temp_dir) / "b_block"
    if temp_path.exists():
        shutil.rmtree(temp_path)
    temp_path.mkdir(parents=True, exist_ok=True)
    export_dir = temp_path / "B_rows.tsv.bgz"

    covariates = hl.literal(covariate_matrix.tolist())
    mt_with_covariates = mt.add_col_index("__ld_col_index")
    mt_with_covariates = mt_with_covariates.annotate_cols(
        __ld_covariates=covariates[hl.int32(mt_with_covariates.__ld_col_index)]
    )

    n_covariates = covariate_matrix.shape[1]
    b_fields = {
        f"b_{i:03d}": hl.format(
            "%.17g",
            hl.agg.sum(
                mt_with_covariates.GT_dosage * mt_with_covariates.__ld_covariates[i]
            ),
        )
        for i in range(n_covariates)
    }
    b_table = mt_with_covariates.annotate_rows(**b_fields).rows()
    b_table = b_table.add_index("__b_row_idx")
    export_columns = ["__b_row_idx", *b_fields.keys()]
    b_table.key_by().select(*export_columns).export(
        str(export_dir),
        header=True,
        parallel="separate_header",
    )

    part_files = [p for p in export_dir.glob("part-*") if p.is_file()]
    if not part_files:
        raise RuntimeError(
            f"B-block export produced no part files under {export_dir}"
        )

    chunks: list[np.ndarray] = []
    for part_file in part_files:
        if part_file.stat().st_size == 0:
            continue
        with gzip.open(part_file, "rt", encoding="utf-8") as handle:
            arr = np.loadtxt(
                handle,
                delimiter="\t",
                dtype=np.float64,
                ndmin=2,
            )
        if arr.size == 0:
            continue
        if arr.shape[1] != n_covariates + 1:
            raise ValueError(
                f"B partition {part_file.name} has {arr.shape[1]} columns; "
                f"expected {n_covariates + 1} (row index + {n_covariates} covariates)"
            )
        chunks.append(arr)

    if not chunks:
        raise RuntimeError(
            f"B-block export produced only empty partitions under {export_dir}"
        )

    combined = np.concatenate(chunks, axis=0)
    order = np.argsort(combined[:, 0].astype(np.int64), kind="stable")
    b_matrix = combined[order, 1:]
    shutil.rmtree(temp_path)
    return np.ascontiguousarray(b_matrix, dtype=np.float64)


DEFAULT_LD_RADIUS_BP = 1_000_000
DEFAULT_A_BLOCK_SIZE = 4096
DEFAULT_A_CHUNK_ROWS = 50_000
DEFAULT_A_MAX_DENSE_GB = 1.0


def _dense_product_gb(n_rows: int, n_cols: int) -> float:
    """Estimate a dense float64 product block size in GiB."""
    return (n_rows * n_cols * 8) / (1024**3)


def _bounded_row_stop(
    row_start: int,
    requested_row_stop: int,
    stops: np.ndarray,
    max_dense_gb: float,
) -> int:
    """Find the largest row stop that stays within the dense-product cap."""
    best_stop: int | None = None
    low = row_start + 1
    high = requested_row_stop

    while low <= high:
        mid = (low + high) // 2
        col_stop = int(stops[row_start:mid].max())
        estimate_gb = _dense_product_gb(mid - row_start, col_stop - row_start)
        if estimate_gb <= max_dense_gb:
            best_stop = mid
            low = mid + 1
        else:
            high = mid - 1

    if best_stop is None:
        col_stop = int(stops[row_start : row_start + 1].max())
        estimate_gb = _dense_product_gb(1, col_stop - row_start)
        raise MemoryError(
            "A-block chunk would exceed the dense intermediate cap even for "
            f"one row: row_start={row_start}, n_cols={col_stop - row_start}, "
            f"estimated_dense_gb={estimate_gb:.3f}, "
            f"max_dense_gb={max_dense_gb:.3f}"
        )

    return best_stop


def compute_a_block_banded(
    mt: hl.MatrixTable,
    variant_index: "VariantIndex",
    out_dir: str | Path,
    radius_bp: int = DEFAULT_LD_RADIUS_BP,
    block_size: int = DEFAULT_A_BLOCK_SIZE,
    chunk_rows: int = DEFAULT_A_CHUNK_ROWS,
    max_dense_gb: float = DEFAULT_A_MAX_DENSE_GB,
) -> dict:
    """Compute and write chunked upper-triangular windowed ``A_k = X X^T``.

    For each chromosome present in ``variant_index``, the MatrixTable is
    filtered to that chromosome, ``X`` is built as a Hail BlockMatrix from
    the imputed ``GT_dosage`` entries, and each row chunk is multiplied only
    by the chromosome column interval needed to cover that chunk's 1 Mb
    upper-triangular windows. Each chunk is then sparsified to row-specific
    intervals and written to ``{out_dir}/chr<C>/chunk_<N>/`` as a Hail
    BlockMatrix directory. The dense chromosome-wide ``p × p`` matrix is
    never materialised.

    Returns a manifest dict describing the written artefacts.
    """
    if chunk_rows <= 0:
        raise ValueError(f"chunk_rows must be positive; got {chunk_rows}")
    if max_dense_gb <= 0:
        raise ValueError(f"max_dense_gb must be positive; got {max_dense_gb}")

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
        stops = j_max_inclusive + 1

        chr_dir = out_path / f"chr{chrom}"
        if chr_dir.exists():
            shutil.rmtree(chr_dir)
        chr_dir.mkdir(parents=True, exist_ok=True)

        mt_chr = mt.filter_rows(mt.locus.contig == chrom)
        x_h: BlockMatrix | None = None

        chunks = []
        chunk_index = 0
        row_start = 0
        while row_start < n_variants:
            requested_row_stop = min(row_start + chunk_rows, n_variants)
            row_stop = _bounded_row_stop(
                row_start,
                requested_row_stop,
                stops,
                max_dense_gb,
            )
            row_indices = list(range(row_start, row_stop))
            chunk_name = f"chunk_{chunk_index:06d}"
            chunk_dir = chr_dir / chunk_name

            col_start = row_start
            col_stop = int(stops[row_start:row_stop].max())
            col_indices = list(range(col_start, col_stop))
            estimated_dense_gb = _dense_product_gb(
                row_stop - row_start,
                col_stop - col_start,
            )

            print(
                "Writing A block "
                f"chr{chrom} {chunk_name}: rows {row_start}-{row_stop}, "
                f"cols {col_start}-{col_stop}, "
                f"estimated dense intermediate {estimated_dense_gb:.3f} GiB",
                flush=True,
            )

            if x_h is None:
                x_h = BlockMatrix.from_entry_expr(
                    mt_chr.GT_dosage,
                    block_size=block_size,
                )
            x_chunk = x_h.filter_rows(row_indices)
            x_window = x_h.filter_rows(col_indices)
            a_chunk = x_chunk @ x_window.T
            a_window = a_chunk.sparsify_row_intervals(
                starts=idx[row_start:row_stop] - col_start,
                stops=stops[row_start:row_stop] - col_start,
                blocks_only=False,
            )
            a_window.write(str(chunk_dir), overwrite=True)

            chunks.append(
                {
                    "name": chunk_name,
                    "directory": f"{out_path.name}/chr{chrom}/{chunk_name}",
                    "row_start": int(row_start),
                    "row_stop": int(row_stop),
                    "n_rows": int(row_stop - row_start),
                    "column_start": int(col_start),
                    "column_stop": int(col_stop),
                    "n_cols": int(col_stop - col_start),
                    "estimated_dense_gb": estimated_dense_gb,
                    "row_index_base": "chromosome",
                    "column_index_base": "chromosome",
                }
            )

            row_start = row_stop
            chunk_index += 1

        chr_manifest = {
            "chromosome": chrom,
            "format": "hail-blockmatrix-row-interval-chunks",
            "radius_bp": radius_bp,
            "block_size": block_size,
            "chunk_rows": chunk_rows,
            "max_dense_gb": max_dense_gb,
            "chunking": "memory-capped up to chunk_rows",
            "n_variants": int(n_variants),
            "n_chunks": len(chunks),
            "sparsification": "row_intervals",
            "window_definition": (
                "upper triangle, same chromosome, pos_j <= pos_i + radius_bp"
            ),
            "max_idx_distance_in_window": max_idx_distance,
            "chunks": chunks,
        }
        (chr_dir / "manifest.json").write_text(
            json.dumps(chr_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        chromosomes_meta[chrom] = {
            "directory": f"{out_path.name}/chr{chrom}",
            "n_variants": int(n_variants),
            "n_chunks": len(chunks),
            "sparsification": "row_intervals",
            "window_definition": (
                "upper triangle, same chromosome, pos_j <= pos_i + radius_bp"
            ),
            "max_idx_distance_in_window": max_idx_distance,
            "block_size": block_size,
            "chunk_rows": chunk_rows,
            "max_dense_gb": max_dense_gb,
            "chunking": "memory-capped up to chunk_rows",
            "chunks": chunks,
        }

    return {
        "directory": out_path.name,
        "format": "hail-blockmatrix-row-interval-chunks",
        "radius_bp": radius_bp,
        "block_size": block_size,
        "chunk_rows": chunk_rows,
        "max_dense_gb": max_dense_gb,
        "chunking": "memory-capped up to chunk_rows",
        "chromosomes": chromosomes_meta,
    }


def read_a_block_chunk(chunk_dir: str | Path) -> np.ndarray:
    """Read a written A_blocks chunk BlockMatrix back as a dense float64 array."""
    return BlockMatrix.read(str(chunk_dir)).to_numpy().astype(np.float64, copy=False)


def write_r_block_chunk(
    dense: np.ndarray,
    starts: np.ndarray,
    stops: np.ndarray,
    out_dir: str | Path,
    block_size: int = DEFAULT_A_BLOCK_SIZE,
) -> None:
    """Write a dense R block as a row-interval BlockMatrix directory.

    Hail BlockMatrix is float64-only (``from_numpy`` upcasts and the
    write/read round-trip always yields float64), so float32 storage is not
    achievable here; the block is persisted as float64.
    `starts`/`stops` are per-row chunk-local column intervals (same convention
    as compute_a_block_banded's sparsify_row_intervals call).
    """
    bm = BlockMatrix.from_numpy(dense.astype(np.float64), block_size=block_size)
    bm = bm.sparsify_row_intervals(
        starts=[int(s) for s in starts],
        stops=[int(s) for s in stops],
        blocks_only=False,
    )
    bm.write(str(out_dir), overwrite=True)
