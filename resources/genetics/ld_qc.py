from pathlib import Path
from typing import TypedDict


type CovariateRow = dict[str, str]
type CovariateTable = dict[str, CovariateRow]


class SampleAlignment(TypedDict):
    covariate_header: list[str]
    sample_counts: dict[str, int]
    final_samples: list[str]


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
DEFAULT_SCHEMA_COLUMNS: list[str] = [
    "Age_numeric",
    "Sex_factor",
]
MISSING_VALUES: set[str] = {"", "NA", "NaN", "nan", "NAN", "N/A", "."}
AUTOSOMES: set[str] = {str(c) for c in range(1, 23)}
VALID_BASES: set[str] = {"A", "C", "G", "T"}


def require_file(path: str | Path, description: str) -> Path:
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing {description}: {path}")
    return path


def read_header(path: str | Path) -> list[str]:
    with Path(path).open("r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()
    if not first_line:
        raise ValueError(f"File has no header: {path}")
    return first_line.split()


def validate_covariate_header(path: str | Path) -> list[str]:
    header = read_header(path)
    missing = [column for column in DEFAULT_SCHEMA_COLUMNS if column not in header]
    if missing:
        raise ValueError(
            "Covariates file is missing required columns for section 15: "
            + ", ".join(missing)
        )
    return header


def read_fam_samples(path: str | Path) -> list[str]:
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


def has_missing_required_covariates(row: CovariateRow) -> bool:
    return any(row[column] in MISSING_VALUES for column in DEFAULT_SCHEMA_COLUMNS)


def build_sample_alignment(
    fam_path: str | Path, covariates_path: str | Path
) -> SampleAlignment:
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

    fam_sample_set = set(fam_samples)
    covariate_sample_set = set(covariates)
    samples_in_both = [iid for iid in fam_samples if iid in covariates]
    final_samples = [
        iid
        for iid in samples_in_both
        if not has_missing_required_covariates(covariates[iid])
    ]

    counts: dict[str, int] = {
        "fam_sample_count": len(fam_samples),
        "covariate_sample_count": len(covariates),
        "samples_in_both_count": len(samples_in_both),
        "fam_without_covariates_count": len(fam_sample_set - covariate_sample_set),
        "covariates_without_fam_count": len(covariate_sample_set - fam_sample_set),
        "missing_required_covariates_count": len(samples_in_both) - len(final_samples),
        "final_sample_count": len(final_samples),
    }

    if not final_samples:
        raise ValueError(
            "No samples remain after FAM/covariate alignment and required covariate filtering"
        )

    return {
        "covariate_header": covariate_header,
        "sample_counts": counts,
        "final_samples": final_samples,
    }


def build_variant_index(bim_path: str | Path) -> VariantIndex:
    bim_path = Path(bim_path)
    seen_keys: set[str] = set()
    duplicates: set[str] = set()
    kept: list[VariantRecord] = []
    counts: dict[str, int] = {
        "total_rows": 0,
        "kept_count": 0,
        "excluded_non_autosomal": 0,
        "excluded_non_biallelic_snp": 0,
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

    kept.sort(key=lambda v: (int(v["chr"]), v["pos"], v["ref"], v["alt"]))
    counts["kept_count"] = len(kept)

    if not kept:
        raise ValueError(
            f"No autosomal biallelic SNPs remain in BIM file: {bim_path}"
        )

    return {"variants": kept, "counts": counts}
