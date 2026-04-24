#!/usr/bin/env python

import argparse
from datetime import datetime, timezone
from pathlib import Path

from ld_io import ensure_dir, write_gzip_text, write_json, write_text
from ld_qc import DEFAULT_SCHEMA_ID, require_file, validate_covariate_header, validate_pcs_header


def parse_args():
    parser = argparse.ArgumentParser(description="Prepare section-15 LD scaffold outputs")
    parser.add_argument("--study-name", required=True)
    parser.add_argument("--bfile", required=True)
    parser.add_argument("--covariates", required=True)
    parser.add_argument("--pcs", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--log-file", required=True)
    return parser.parse_args()


def main():
    args = parse_args()

    output_dir = ensure_dir(args.output_dir)
    blocks_dir = ensure_dir(output_dir / "blocks")

    for suffix in (".bed", ".bim", ".fam"):
        require_file(f"{args.bfile}{suffix}", f"cleaned genotype input {suffix}")

    covariate_columns = validate_covariate_header(require_file(args.covariates, "covariates input"))
    pc_columns = validate_pcs_header(require_file(args.pcs, "genotype PCs input"))

    manifest = {
        "study_name": args.study_name,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "module": "15a",
        "status": "scaffold",
        "genome_build": "GRCh37",
        "autosomes_only": True,
        "input_bfile": args.bfile,
        "covariate_schema": {
            "schema_id": DEFAULT_SCHEMA_ID,
            "required_columns": ["Age_numeric", "Sex_factor", "PC1-PC10"],
            "covariates_file": args.covariates,
            "pcs_file": args.pcs,
        },
        "notes": [
            "This scaffold validates required section-15 inputs and writes placeholder outputs.",
            "Production A/B/D matrix generation is not implemented yet.",
        ],
    }

    variant_stub = "chr\tpos\tref\talt\tvariant_id\n"
    diagnostics = [
        "Section 15 cohort scaffold created successfully.",
        f"Covariate columns detected: {', '.join(covariate_columns)}",
        f"PC header columns detected: {', '.join(pc_columns[:12])}",
        "Next implementation step: replace placeholder outputs with Hail/BlockMatrix cohort statistics.",
    ]

    write_json(output_dir / "manifest.json", manifest)
    write_gzip_text(output_dir / "variants.tsv.gz", variant_stub)
    write_text(output_dir / "qc_report.txt", "\n".join(diagnostics) + "\n")
    write_text(blocks_dir / ".gitkeep", "")


if __name__ == "__main__":
    main()
