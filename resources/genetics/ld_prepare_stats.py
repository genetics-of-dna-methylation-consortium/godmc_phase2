#!/usr/bin/env python

import argparse
from datetime import datetime, timezone
from typing import Iterable

from ld_io import JsonValue, ensure_dir, write_gzip_lines, write_json, write_text
from ld_qc import (
    DEFAULT_SCHEMA_ID,
    VariantRecord,
    build_sample_alignment,
    build_variant_index,
    require_file,
    validate_covariate_header,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Prepare section-15 LD scaffold outputs")
    parser.add_argument("--study-name", required=True)
    parser.add_argument("--bfile", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    return parser.parse_args()


VARIANT_TSV_HEADER = "chr\tpos\tref\talt\tvariant_id"


def _variant_lines(variants: list[VariantRecord]) -> Iterable[str]:
    yield VARIANT_TSV_HEADER
    for variant in variants:
        yield (
            f"{variant['chr']}\t{variant['pos']}\t{variant['ref']}\t"
            f"{variant['alt']}\t{variant['variant_id']}"
        )


def main() -> None:
    args = parse_args()

    output_dir = ensure_dir(args.output_dir)
    blocks_dir = ensure_dir(output_dir / "blocks")

    for suffix in (".bed", ".bim", ".fam"):
        require_file(f"{args.bfile}{suffix}", f"cleaned genotype input {suffix}")

    covariate_columns = validate_covariate_header(
        require_file(args.covariates, "covariates input")
    )
    sample_alignment = build_sample_alignment(f"{args.bfile}.fam", args.covariates)
    sample_counts = sample_alignment["sample_counts"]

    variant_index = build_variant_index(f"{args.bfile}.bim")
    variant_counts = variant_index["counts"]

    manifest: dict[str, JsonValue] = {
        "study_name": args.study_name,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "module": "15a",
        "status": "scaffold",
        "genome_build": "GRCh37",
        "autosomes_only": True,
        "input_bfile": args.bfile,
        "sample_alignment": {
            "sample_order": "cleaned section-02 FAM order after required covariate filtering",
            "counts": sample_counts,
        },
        "variant_index": {
            "schema_version": "v0.2-index-only",
            "columns": ["chr", "pos", "ref", "alt", "variant_id"],
            "deferred_columns": ["n_nonmissing", "n_imputed", "genotype_mean"],
            "filters_applied": [
                "autosomes 1-22",
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
            "covariates_file": args.covariates,
        },
        "notes": [
            "This scaffold validates required section-15 inputs and writes placeholder outputs.",
            "The default first-release covariate schema excludes cohort-specific genotype PCs.",
            "Production A/B/D matrix generation is not implemented yet.",
            "variants.tsv.gz currently contains index columns only; per-variant genotype "
            "diagnostics will be added when .bed reading is implemented.",
        ],
    }

    diagnostics: list[str] = [
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
        "Samples dropped for missing required covariates: "
        f"{sample_counts['missing_required_covariates_count']}",
        f"Final section-15 sample count: {sample_counts['final_sample_count']}",
        f"BIM rows scanned: {variant_counts['total_rows']}",
        f"Variants kept (autosomal biallelic SNPs): {variant_counts['kept_count']}",
        "Variants excluded as non-autosomal: "
        f"{variant_counts['excluded_non_autosomal']}",
        "Variants excluded as non-biallelic SNP: "
        f"{variant_counts['excluded_non_biallelic_snp']}",
        "Next implementation step: read genotype calls from .bed for D_k export "
        "and per-variant diagnostics.",
    ]

    write_json(output_dir / "manifest.json", manifest)
    write_gzip_lines(output_dir / "variants.tsv.gz", _variant_lines(variant_index["variants"]))
    write_text(output_dir / "qc_report.txt", "\n".join(diagnostics) + "\n")
    write_text(blocks_dir / ".gitkeep", "")


if __name__ == "__main__":
    main()
