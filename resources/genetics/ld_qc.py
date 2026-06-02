from pathlib import Path
from typing import TypedDict

import numpy as np


type CovariateRow = dict[str, str]
type CovariateTable = dict[str, CovariateRow]


class SampleAlignment(TypedDict):
    covariate_header: list[str]
    sample_counts: dict[str, int]
    final_samples: list[str]
    covariates: CovariateTable


class VariantRecord(TypedDict):
    chr: str
    pos: int
    ref: str
    alt: str
    variant_id: str


class VariantIndex(TypedDict):
    variants: list[VariantRecord]
    counts: dict[str, int]


DEFAULT_SCHEMA_ID = "intercept_age_sex"
DEFAULT_SCHEMA_COLUMNS = [
    "Age_numeric",
    "Sex_factor",
]
COVARIATE_MATRIX_COLUMNS = ["intercept", "Age_numeric", "Sex_factor"]
SEX_FACTOR_RECODE = {"M": 1.0, "F": 2.0}
AUTOSOMES = {str(c) for c in range(1, 23)}
VALID_BASES = {"A", "C", "G", "T"}
MHC_CHROMOSOME = "6"
MHC_START_BP = 28477797
MHC_END_BP = 33448354


def require_file(path: str | Path, description: str) -> Path:
    """Return an existing file path or raise a clear missing-input error."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing {description}: {path}")
    return path


def read_fam_samples(path: str | Path) -> list[str]:
    """Read sample IIDs from the second column of a PLINK FAM file."""
    samples: list[str] = []
    seen: set[str] = set()
    duplicates: set[str] = set()
    with Path(path).open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            row = line.split()
            if not row:
                continue
            if len(row) < 2:
                raise ValueError(
                    f"FAM row {line_number} has fewer than two columns: {path}"
                )
            iid = row[1]
            if iid in seen:
                duplicates.add(iid)
            seen.add(iid)
            samples.append(iid)
    if duplicates:
        raise ValueError(f"Duplicate IID values in FAM file: {len(duplicates)}")
    return samples


def read_covariates(path: str | Path) -> tuple[list[str], CovariateTable]:
    """Read the whitespace-delimited covariate table keyed by IID."""
    path = Path(path)
    with path.open("r", encoding="utf-8") as handle:
        header_line = handle.readline().strip()
        if not header_line:
            raise ValueError(f"File has no header: {path}")
        header = header_line.split()
        if "IID" not in header:
            raise ValueError(f"Covariates file is missing required IID column: {path}")

        iid_index = header.index("IID")
        covariates: CovariateTable = {}
        duplicates: set[str] = set()
        for line_number, line in enumerate(handle, start=2):
            row = line.split()
            if not row:
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Covariates row {line_number} has {len(row)} columns; "
                    f"expected {len(header)}: {path}"
                )
            iid = row[iid_index]
            if iid in covariates:
                duplicates.add(iid)
            covariates[iid] = dict(zip(header, row))

    if duplicates:
        raise ValueError(f"Duplicate IID values in covariates file: {len(duplicates)}")
    return header, covariates


def build_sample_alignment(
    fam_path: str | Path, covariates_path: str | Path
) -> SampleAlignment:
    """Align FAM samples to covariate rows in cleaned genotype sample order."""
    fam_samples = read_fam_samples(fam_path)
    covariate_header, covariates = read_covariates(covariates_path)
    missing_columns = [
        column for column in DEFAULT_SCHEMA_COLUMNS if column not in covariate_header
    ]
    if missing_columns:
        raise ValueError(
            "Covariates file is missing required columns for section 15: "
            + ", ".join(missing_columns)
        )

    # Assemble list of individuals in both FAM and covariate files
    fam_sample_set = set(fam_samples)
    covariate_sample_set = set(covariates)
    samples_in_both = [iid for iid in fam_samples if iid in covariates]
    final_samples = samples_in_both

    counts = {
        "fam_sample_count": len(fam_samples),
        "covariate_sample_count": len(covariates),
        "samples_in_both_count": len(samples_in_both),
        "fam_without_covariates_count": len(fam_sample_set - covariate_sample_set),
        "covariates_without_fam_count": len(covariate_sample_set - fam_sample_set),
        "final_sample_count": len(final_samples),
    }

    if not final_samples:
        raise ValueError(
            "No samples remain after FAM/covariate alignment"
        )

    return {
        "covariate_header": covariate_header,
        "sample_counts": counts,
        "final_samples": final_samples,
        "covariates": covariates,
    }


def build_covariate_matrix(
    final_samples: list[str], covariates: CovariateTable
) -> np.ndarray:
    """Build the frozen section-15 covariate matrix for aligned samples."""
    if not final_samples:
        raise ValueError("Cannot build covariate matrix from empty sample list")
    matrix = np.empty((len(final_samples), len(COVARIATE_MATRIX_COLUMNS)), dtype=np.float64)
    for i, iid in enumerate(final_samples):
        if iid not in covariates:
            raise KeyError(f"Sample {iid} missing from covariate table")
        row = covariates[iid]
        try:
            age = float(row["Age_numeric"])
        except (KeyError, ValueError) as exc:
            raise ValueError(
                f"Sample {iid} has invalid Age_numeric "
                f"'{row.get('Age_numeric')}'"
            ) from exc
        sex_raw = row.get("Sex_factor", "")
        if sex_raw not in SEX_FACTOR_RECODE:
            raise ValueError(
                f"Sample {iid} has unrecognised Sex_factor '{sex_raw}'; "
                f"expected one of {sorted(SEX_FACTOR_RECODE)}"
            )
        matrix[i, 0] = 1.0
        matrix[i, 1] = age
        matrix[i, 2] = SEX_FACTOR_RECODE[sex_raw]
    return matrix


def build_variant_index(
    bim_path: str | Path, chromosome: str | None = None
) -> VariantIndex:
    """Build canonical autosomal variant records from a cleaned PLINK BIM file."""
    bim_path = Path(bim_path)
    if chromosome is not None and chromosome not in AUTOSOMES:
        raise ValueError(
            f"Invalid chromosome filter '{chromosome}'; "
            f"expected one of {sorted(AUTOSOMES, key=int)} or None for all autosomes"
        )
    seen_keys: set[str] = set()
    duplicates: set[str] = set()
    kept: list[VariantRecord] = []
    counts: dict[str, int] = {
        "total_rows": 0,
        "kept_count": 0,
        "excluded_non_autosomal": 0,
        "excluded_other_chromosome": 0,
        "excluded_non_biallelic_snp": 0,
        "excluded_mhc_region": 0,
    }

    with bim_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            row = line.split()
            if not row:
                continue
            if len(row) < 6:
                raise ValueError(
                    f"BIM row {line_number} has fewer than six columns: {bim_path}"
                )
            counts["total_rows"] += 1
            chr_raw, _vid, _cm, pos_raw, allele1, allele2 = row[:6]

            if chr_raw not in AUTOSOMES:
                counts["excluded_non_autosomal"] += 1
                continue

            if chromosome is not None and chr_raw != chromosome:
                counts["excluded_other_chromosome"] += 1
                continue

            ref = allele2.upper()
            alt = allele1.upper()
            if ref not in VALID_BASES or alt not in VALID_BASES or ref == alt:
                counts["excluded_non_biallelic_snp"] += 1
                continue

            try:
                pos = int(pos_raw)
            except ValueError as exc:
                raise ValueError(
                    f"BIM row {line_number} has non-integer position '{pos_raw}': "
                    f"{bim_path}"
                ) from exc
            if pos <= 0:
                raise ValueError(
                    f"BIM row {line_number} has non-positive position {pos}: {bim_path}"
                )

            if (
                chr_raw == MHC_CHROMOSOME
                and MHC_START_BP <= pos <= MHC_END_BP
            ):
                counts["excluded_mhc_region"] += 1
                continue

            variant_id = f"{chr_raw}:{pos}:{ref}:{alt}"
            if variant_id in seen_keys:
                duplicates.add(variant_id)
                continue
            seen_keys.add(variant_id)
            kept.append(
                {
                    "chr": chr_raw,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "variant_id": variant_id,
                }
            )

    if duplicates:
        raise ValueError(
            f"Duplicate canonical variant keys in BIM file: {len(duplicates)} "
            f"({bim_path})"
        )

    sorted_kept = sorted(
        kept, key=lambda v: (int(v["chr"]), v["pos"], v["ref"], v["alt"])
    )
    if kept != sorted_kept:
        first_diff = next(
            i for i, (observed, expected) in enumerate(zip(kept, sorted_kept))
            if observed != expected
        )
        raise ValueError(
            "BIM variants must already be sorted by chr, position, ref, alt for "
            "section 15 Hail row alignment; first mismatch at kept variant row "
            f"{first_diff}: observed {kept[first_diff]['variant_id']}, "
            f"expected {sorted_kept[first_diff]['variant_id']}"
        )

    kept = sorted_kept
    counts["kept_count"] = len(kept)

    if not kept:
        scope = (
            f"chromosome {chromosome}"
            if chromosome is not None
            else "autosomes 1-22"
        )
        raise ValueError(
            f"No biallelic SNPs remain on {scope} in BIM file: {bim_path}"
        )

    return {"variants": kept, "counts": counts}
