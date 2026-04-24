from pathlib import Path


DEFAULT_SCHEMA_ID = "intercept_age_sex_pcs1_10"
DEFAULT_SCHEMA_COLUMNS = [
    "Age_numeric",
    "Sex_factor",
]


def require_file(path, description):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"Missing {description}: {path}")
    return path


def read_header(path):
    with Path(path).open("r", encoding="utf-8") as handle:
        first_line = handle.readline().strip()
    if not first_line:
        raise ValueError(f"File has no header: {path}")
    return first_line.split()


def validate_covariate_header(path):
    header = read_header(path)
    missing = [column for column in DEFAULT_SCHEMA_COLUMNS if column not in header]
    if missing:
        raise ValueError(
            "Covariates file is missing required columns for section 15: "
            + ", ".join(missing)
        )
    return header


def validate_pcs_header(path):
    row = read_header(path)
    if len(row) < 12:
        raise ValueError(
            "PC file does not contain IID plus at least 10 PCs required for section 15"
        )
    return row
