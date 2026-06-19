#!/usr/bin/env python

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ld_hail import (
    DEFAULT_A_BLOCK_SIZE,
    DEFAULT_A_CHUNK_ROWS,
    DEFAULT_A_MAX_DENSE_GB,
    DEFAULT_N_PARTITIONS,
    compute_a_block_banded,
    compute_b_block,
    init_hail,
    load_genotype_matrixtable,
    prepare_for_cross_products,
)
from ld_qc import (
    COVARIATE_MATRIX_COLUMNS,
    DEFAULT_SCHEMA_ID,
    build_covariate_matrix,
    build_sample_alignment,
    build_variant_index,
    require_file,
)


CHROMOSOME_FILTER_ALL = "all"


def parse_args() -> argparse.Namespace:
    """Parse command-line options for cohort-side LD scaffold preparation."""
    parser = argparse.ArgumentParser(description="Prepare section-15 LD scaffold outputs")
    parser.add_argument("--study-name", required=True)
    parser.add_argument("--bfile", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument(
        "--hail-partitions",
        type=int,
        default=DEFAULT_N_PARTITIONS,
        help=f"Number of Hail import partitions (default: {DEFAULT_N_PARTITIONS}).",
    )
    parser.add_argument(
        "--hail-local-cores",
        type=int,
        default=None,
        help="Optional local Spark core count for Hail, e.g. SLURM_CPUS_ON_NODE.",
    )
    parser.add_argument(
        "--hail-driver-memory-gb",
        type=int,
        default=None,
        help="Optional Spark driver memory in GB for Hail.",
    )
    parser.add_argument(
        "--hail-tmp-dir",
        default=None,
        help="Optional Hail temporary directory (default: output-dir/hail_tmp).",
    )
    parser.add_argument(
        "--a-block-size",
        type=int,
        default=DEFAULT_A_BLOCK_SIZE,
        help=f"Hail BlockMatrix block size for A chunks (default: {DEFAULT_A_BLOCK_SIZE}).",
    )
    parser.add_argument(
        "--a-chunk-rows",
        type=int,
        default=DEFAULT_A_CHUNK_ROWS,
        help=f"Number of variant rows per A chunk (default: {DEFAULT_A_CHUNK_ROWS}).",
    )
    parser.add_argument(
        "--a-max-dense-gb",
        type=float,
        default=DEFAULT_A_MAX_DENSE_GB,
        help=(
            "Maximum estimated dense intermediate size per A chunk in GiB "
            f"(default: {DEFAULT_A_MAX_DENSE_GB})."
        ),
    )
    parser.add_argument(
        "--chromosome",
        default=CHROMOSOME_FILTER_ALL,
        help=(
            f"Restrict to a single autosome (1-22), or '{CHROMOSOME_FILTER_ALL}' for "
            f"all autosomes (default: {CHROMOSOME_FILTER_ALL})."
        ),
    )
    return parser.parse_args()


def resolve_chromosome_filter(value: str) -> str | None:
    """Convert the CLI chromosome value into an optional autosome filter."""
    if value == CHROMOSOME_FILTER_ALL:
        return None
    return value


def _export_variants_tsv_gz(path: Path, diagnostics_table) -> None:
    """Write canonical variant index rows + genotype diagnostics via Hail."""
    diagnostics_table.key_by().select(
        "chr",
        "pos",
        "ref",
        "alt",
        "variant_id",
        "n_nonmissing",
        "n_imputed",
        "genotype_mean",
    ).export(str(path), header=True)


def log_step(message: str) -> None:
    """Emit a flushed progress marker into the section log."""
    print(f"[section15a] {message}", flush=True)


def main() -> None:
    """Create section-15 cohort scaffold outputs from cleaned pipeline inputs."""
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    blocks_dir = output_dir / "blocks"
    blocks_dir.mkdir(parents=True, exist_ok=True)
    hail_tmp_dir = (
        Path(args.hail_tmp_dir) if args.hail_tmp_dir else output_dir / "hail_tmp"
    )

    log_step("Checking section-02 genotype and covariate inputs")
    for suffix in (".bed", ".bim", ".fam"):
        require_file(f"{args.bfile}{suffix}", f"cleaned genotype input {suffix}")

    require_file(args.covariates, "covariates input")
    log_step("Building sample alignment and canonical variant index")
    sample_alignment = build_sample_alignment(f"{args.bfile}.fam", args.covariates)
    sample_counts = sample_alignment["sample_counts"]
    covariate_columns = sample_alignment["covariate_header"]

    chromosome_filter = resolve_chromosome_filter(args.chromosome)
    variant_index = build_variant_index(
        f"{args.bfile}.bim", chromosome=chromosome_filter
    )
    variant_counts = variant_index["counts"]

    covariate_matrix = build_covariate_matrix(
        sample_alignment["final_samples"], sample_alignment["covariates"]
    )
    d_matrix = covariate_matrix.T @ covariate_matrix
    np.save(output_dir / "D.npy", d_matrix, allow_pickle=False)
    d_rank = int(np.linalg.matrix_rank(d_matrix))
    d_condition_number = float(np.linalg.cond(d_matrix))

    hail_log = Path(args.log_file).parent / "hail.log"
    log_step("Initialising Hail")
    hail_meta = init_hail(
        hail_log,
        n_partitions=args.hail_partitions,
        local_cores=args.hail_local_cores,
        driver_memory_gb=args.hail_driver_memory_gb,
        tmp_dir=hail_tmp_dir,
    )
    log_step("Loading cleaned genotype MatrixTable")
    mt = load_genotype_matrixtable(
        bfile=args.bfile,
        chromosome=chromosome_filter,
        final_samples=sample_alignment["final_samples"],
        variant_index=variant_index,
        n_partitions=args.hail_partitions,
    )
    log_step("Computing genotype diagnostics and mean-imputed dosage entries")
    mt_imputed, genotype_diagnostics = prepare_for_cross_products(
        mt,
        n_samples=sample_counts["final_sample_count"],
    )

    log_step("Computing B = X^T C using row aggregations")
    b_matrix = compute_b_block(mt_imputed, covariate_matrix, temp_dir=hail_tmp_dir)
    np.save(output_dir / "B.npy", b_matrix, allow_pickle=False)
    b_frobenius = float(np.linalg.norm(b_matrix))
    log_step("Exporting variant diagnostics before A-block computation")
    _export_variants_tsv_gz(output_dir / "variants.tsv.gz", genotype_diagnostics)

    log_step("Computing memory-capped A = X X^T row-interval chunks")
    a_blocks_meta = compute_a_block_banded(
        mt_imputed,
        variant_index,
        out_dir=output_dir / "A_blocks",
        block_size=args.a_block_size,
        chunk_rows=args.a_chunk_rows,
        max_dense_gb=args.a_max_dense_gb,
    )

    manifest = {
        "study_name": args.study_name,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "module": "15a",
        "status": "scaffold",
        "genome_build": "GRCh37",
        "autosomes_only": True,
        "input_bfile": args.bfile,
        "sample_alignment": {
            "sample_order": "cleaned section-02 FAM order after covariate alignment",
            "counts": sample_counts,
        },
        "variant_index": {
            "schema_version": "v0.3-with-genotype-stats",
            "columns": [
                "chr",
                "pos",
                "ref",
                "alt",
                "variant_id",
                "n_nonmissing",
                "n_imputed",
                "genotype_mean",
            ],
            "chromosome_filter": args.chromosome,
            "filters_applied": [
                f"chromosome {chromosome_filter}"
                if chromosome_filter is not None
                else "autosomes 1-22",
                "biallelic SNPs (ref/alt in {A,C,G,T}, ref != alt)",
                "MHC region chr6:28477797-33448354 excluded (GRCh37/hg19)",
            ],
            "ref_alt_convention": "ref = .bim column 6 (A2), alt = .bim column 5 (A1)",
            "duplicate_key_policy": "hard fail on duplicate chr:pos:ref:alt",
            "sort_order": "chr (numeric), pos, ref, alt",
            "counts": variant_counts,
            "n_samples_used": sample_counts["final_sample_count"],
        },
        "covariate_schema": {
            "schema_id": DEFAULT_SCHEMA_ID,
            "required_columns": [],
            "matrix_columns": COVARIATE_MATRIX_COLUMNS,
            "sex_factor_recode": {},
            "covariates_file": args.covariates,
        },
        "hail": hail_meta,
        "B_block": {
            "filename": "B.npy",
            "shape": list(b_matrix.shape),
            "dtype": str(b_matrix.dtype),
            "rows": "variants.tsv.gz row order",
            "columns": COVARIATE_MATRIX_COLUMNS,
            "format": "numpy-npy-dense",
        },
        "A_blocks": a_blocks_meta,
        "notes": [
            "D.npy is the dense covariate cross-product C^T C in float64.",
            "B.npy is the dense X^T C cross-product (variants x covariates), float64.",
            "A_blocks/chr<C>/chunk_<N>/ contains per-chromosome chunked "
            "upper-triangular row-interval X X^T as Hail BlockMatrix directories "
            "within radius_bp physical distance.",
            "The covariate schema is intercept-only (grand-mean centring): it "
            "matches the GoDMC mQTL estimand, where every covariate is regressed "
            "out of methylation in section 03, never out of the genotypes.",
            "Per-variant genotype diagnostics are computed via Hail and written to "
            "variants.tsv.gz; missing calls are mean-imputed before downstream cross-products.",
        ],
    }

    qc_lines = [
        "Section 15 cohort scaffold created successfully.",
        f"Covariate columns detected: {', '.join(covariate_columns)}",
        f"FAM samples: {sample_counts['fam_sample_count']}",
        f"Covariate samples: {sample_counts['covariate_sample_count']}",
        "Samples present in both FAM and covariates: "
        f"{sample_counts['samples_in_both_count']}",
        "FAM samples without covariates: "
        f"{sample_counts['fam_without_covariates_count']}",
        "Covariate samples absent from FAM: "
        f"{sample_counts['covariates_without_fam_count']}",
        f"Final section-15 sample count: {sample_counts['final_sample_count']}",
        f"Chromosome filter: {args.chromosome}",
        f"BIM rows scanned: {variant_counts['total_rows']}",
        f"Variants kept (autosomal biallelic SNPs): {variant_counts['kept_count']}",
        "Variants excluded as non-autosomal: "
        f"{variant_counts['excluded_non_autosomal']}",
        "Variants excluded as off-target chromosome: "
        f"{variant_counts['excluded_other_chromosome']}",
        "Variants excluded as non-biallelic SNP: "
        f"{variant_counts['excluded_non_biallelic_snp']}",
        "Variants excluded in MHC region chr6:28477797-33448354: "
        f"{variant_counts['excluded_mhc_region']}",
        f"D matrix shape: {d_matrix.shape[0]} x {d_matrix.shape[1]}",
        f"D matrix rank: {d_rank} (expected {len(COVARIATE_MATRIX_COLUMNS)})",
        f"D matrix condition number: {d_condition_number:.6g}",
        f"Hail version: {hail_meta['hail_version']}",
        "Mean-imputed genotype calls: see variants.tsv.gz n_imputed column",
        f"B matrix shape: {b_matrix.shape[0]} x {b_matrix.shape[1]}",
        f"B Frobenius norm: {b_frobenius:.6g}",
        f"A_blocks radius_bp: {a_blocks_meta['radius_bp']}",
        f"A_blocks block_size: {a_blocks_meta['block_size']}",
        f"A_blocks chunk_rows: {a_blocks_meta['chunk_rows']}",
        f"A_blocks max_dense_gb: {a_blocks_meta['max_dense_gb']}",
        "A_blocks per-chromosome variant counts: "
        + ", ".join(
            f"chr{c}={meta['n_variants']}"
            for c, meta in a_blocks_meta["chromosomes"].items()
        ),
        "A_blocks per-chromosome chunk counts: "
        + ", ".join(
            f"chr{c}={meta['n_chunks']}"
            for c, meta in a_blocks_meta["chromosomes"].items()
        ),
        "Next implementation step: pooled aggregation in 15b.",
    ]

    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "qc_report.txt").write_text(
        "\n".join(qc_lines) + "\n", encoding="utf-8"
    )
    (blocks_dir / ".gitkeep").write_text("", encoding="utf-8")


if __name__ == "__main__":
    main()
