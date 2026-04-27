#!/usr/bin/env python

import argparse
import gzip
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from ld_qc import (
    COVARIATE_MATRIX_COLUMNS,
    DEFAULT_SCHEMA_ID,
    VariantRecord,
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


VARIANT_TSV_HEADER = "chr\tpos\tref\talt\tvariant_id"


def _write_variants_tsv_gz(path: Path, variants: list[VariantRecord]) -> None:
    """Write canonical variant index rows to a gzipped TSV file."""
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(VARIANT_TSV_HEADER + "\n")
        for v in variants:
            handle.write(
                f"{v['chr']}\t{v['pos']}\t{v['ref']}\t"
                f"{v['alt']}\t{v['variant_id']}\n"
            )


def main() -> None:
    """Create section-15 cohort scaffold outputs from cleaned pipeline inputs."""
    args = parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    blocks_dir = output_dir / "blocks"
    blocks_dir.mkdir(parents=True, exist_ok=True)

    for suffix in (".bed", ".bim", ".fam"):
        require_file(f"{args.bfile}{suffix}", f"cleaned genotype input {suffix}")

    require_file(args.covariates, "covariates input")
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
            "schema_version": "v0.2-index-only",
            "columns": ["chr", "pos", "ref", "alt", "variant_id"],
            "deferred_columns": ["n_nonmissing", "n_imputed", "genotype_mean"],
            "chromosome_filter": args.chromosome,
            "filters_applied": [
                f"chromosome {chromosome_filter}"
                if chromosome_filter is not None
                else "autosomes 1-22",
                "biallelic SNPs (ref/alt in {A,C,G,T}, ref != alt)",
            ],
            "ref_alt_convention": "ref = .bim column 6 (A2), alt = .bim column 5 (A1)",
            "duplicate_key_policy": "hard fail on duplicate chr:pos:ref:alt",
            "sort_order": "chr (numeric), pos, ref, alt",
            "counts": variant_counts,
        },
        "covariate_schema": {
            "schema_id": DEFAULT_SCHEMA_ID,
            "required_columns": ["Age_numeric", "Sex_factor"],
            "matrix_columns": COVARIATE_MATRIX_COLUMNS,
            "sex_factor_recode": {"M": 1.0, "F": 2.0},
            "covariates_file": args.covariates,
        },
        "notes": [
            "D.npy is the dense covariate cross-product C^T C in float64.",
            "The default first-release covariate schema excludes cohort-specific genotype PCs.",
            "B/A block exports and per-variant genotype diagnostics are not implemented yet.",
            "variants.tsv.gz currently contains index columns only; per-variant genotype "
            "diagnostics will be added when .bed reading is implemented.",
        ],
    }

    diagnostics = [
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
        f"D matrix shape: {d_matrix.shape[0]} x {d_matrix.shape[1]}",
        f"D matrix rank: {d_rank} (expected {len(COVARIATE_MATRIX_COLUMNS)})",
        f"D matrix condition number: {d_condition_number:.6g}",
        "Next implementation step: read genotype calls from .bed for B_k/A_k export "
        "and per-variant diagnostics.",
    ]

    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    _write_variants_tsv_gz(output_dir / "variants.tsv.gz", variant_index["variants"])
    (output_dir / "qc_report.txt").write_text(
        "\n".join(diagnostics) + "\n", encoding="utf-8"
    )
    (blocks_dir / ".gitkeep").write_text("", encoding="utf-8")


if __name__ == "__main__":
    main()
