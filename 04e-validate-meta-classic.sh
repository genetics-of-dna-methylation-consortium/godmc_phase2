#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

validation_cpg="${positive_control_cpg:-cg07959070}"
validation_out="${section_04_dir}/meta_classic_validation_positive_control"
meta_inputs="${section_04_dir}/meta_inputs"
run_out="${validation_out}/run"
selected_covariates="${validation_out}/selected_covariates.tsv"
ph_id_inc="${validation_out}/positive_control_cpg.txt"
reference_file="${HASE_REF_FILE:-${light_hase}/data/ref-hrc.ref.gz}"
plink_positive_control="${section_03_dir}/positive_control_transformed_${validation_cpg}.PHENO1.glm.linear.gz"

if [ ! -f "${reference_file}" ] && [ -f "${hase}/data/ref-hrc.ref.gz" ]; then
    reference_file="${hase}/data/ref-hrc.ref.gz"
fi

mkdir -p "${validation_out}" "${run_out}"
exec &> >(tee "${validation_out}/log.txt")
print_version

echo "Validating 04e meta inputs with light_hase meta-classic"
echo "Study: ${study_name}"
echo "Positive control CpG: ${validation_cpg}"
echo "Meta inputs: ${meta_inputs}"
echo "Reference file: ${reference_file}"
echo "PLINK positive-control file: ${plink_positive_control}"
echo "Output: ${validation_out}"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

check_file() {
    if [ ! -f "$1" ]; then
        fail "Missing required file: $1"
    fi
}

check_dir() {
    if [ ! -d "$1" ]; then
        fail "Missing required directory: $1"
    fi
}

check_file "${meta_inputs}/part_dev/${study_name}_a_cov.npy"
check_file "${meta_inputs}/part_dev/${study_name}_b_cov.npy"
check_file "${meta_inputs}/part_dev/${study_name}_C.npy"
check_file "${meta_inputs}/part_dev/${study_name}_a_test.npy"
check_file "${meta_inputs}/part_dev/${study_name}_metadata.npy"

check_dir "${meta_inputs}/use_data"
check_dir "${meta_inputs}/use_data/genotype"
check_dir "${meta_inputs}/use_data/individuals"
check_dir "${meta_inputs}/use_data/phenotypes"
check_dir "${meta_inputs}/use_data/probes"
check_dir "${meta_inputs}/mapping"

if ! find "${meta_inputs}/mapping" -maxdepth 1 -type f -name "*.npy" | grep -q .; then
    fail "Missing mapper npy files in: ${meta_inputs}/mapping"
fi

check_file "${reference_file}"
check_file "${plink_positive_control}"

printf "ID\n%s\n" "${validation_cpg}" > "${ph_id_inc}"

echo "Python runner: ${PYTHON_RUNNER[*]}"

printf "%s\t%s_intercept\n" "${study_name}" "${study_name}" > "${selected_covariates}"
echo "Wrote selected covariates: ${selected_covariates}"
echo "Selected covariates: ${study_name} ${study_name}_intercept"

echo "Running light_hase meta-classic"

"${PYTHON_RUNNER[@]}" "${light_hase}/hase.py" \
    -mode meta-classic \
    -study_name "${study_name}" \
    -g "${meta_inputs}/use_data" \
    -ph "${meta_inputs}/use_data/phenotypes" \
    -derivatives "${meta_inputs}/part_dev" \
    -mapper "${meta_inputs}/mapping" \
    -ph_id_inc "${ph_id_inc}" \
    -encoded 1 \
    --selected-covariates "${selected_covariates}" \
    -ref_name ref-hrc \
    -o "${run_out}" \
    -thr 0 \
    -thr_full_log 0 \
    -max-missingness-rate 1 \
    -cluster n

echo "Combining feather outputs and writing gzip-compressed CSV files"

"${PYTHON_RUNNER[@]}" - \
    "${run_out}" \
    "${validation_out}" \
    "${study_name}" \
    "${validation_cpg}" \
    "${reference_file}" <<'PY'
import glob
import os
import sys

import pandas as pd

run_out, validation_out, study_name, cpg, reference_file = sys.argv[1:6]


def read_feathers(pattern, label):
    files = sorted(glob.glob(pattern))
    if not files:
        raise SystemExit("No {} feather files found with pattern: {}".format(label, pattern))
    frames = []
    for path in files:
        frames.append(pd.read_feather(path))
    result = pd.concat(frames, ignore_index=True)
    if result.empty:
        raise SystemExit("{} feather files were found but combined data is empty".format(label))
    return result, files


def pick_column(columns, candidates, required=True):
    for candidate in candidates:
        if candidate in columns:
            return candidate
    if required:
        raise SystemExit(
            "Reference file {} is missing one of these columns: {}".format(
                reference_file, ", ".join(candidates)))
    return None


def load_reference(reference_file):
    compression = "gzip" if reference_file.endswith(".gz") else None
    ref = pd.read_csv(reference_file, delim_whitespace=True, compression=compression)
    ref = ref.reset_index(drop=True)

    id_col = pick_column(ref.columns, ["ID", "id", "variant", "SNP"])
    allele1_col = pick_column(ref.columns, ["str_allele1", "allele1", "A1", "effect_allele"])
    allele2_col = pick_column(ref.columns, ["str_allele2", "allele2", "A2", "non_effect_allele"])
    chr_col = pick_column(ref.columns, ["CHR", "#CHROM", "chromosome", "chr"], required=False)
    bp_col = pick_column(ref.columns, ["bp", "BP", "pos", "position"], required=False)

    return ref, id_col, allele1_col, allele2_col, chr_col, bp_col


def annotate_variants(df, ref, id_col, allele1_col, allele2_col, chr_col, bp_col, label):
    if "variant_index" not in df.columns:
        raise SystemExit("{} results do not contain a variant_index column".format(label))

    variant_index = df["variant_index"].astype("int64")
    if variant_index.min() < 0 or variant_index.max() >= ref.shape[0]:
        raise SystemExit(
            "{} variant_index values are outside reference row range 0-{}".format(
                label, ref.shape[0] - 1))

    annotated = df.copy()
    annotated["ID"] = variant_index.map(ref[id_col])
    # HASE/light_hase decode PLINK .bed genotypes as .bim allele2 dosage.
    # Therefore the HASE beta is relative to str_allele2, while str_allele1
    # is the other allele, even though PLINK often labels A2 as "other".
    annotated["hase_beta_allele"] = variant_index.map(ref[allele2_col])
    annotated["hase_other_allele"] = variant_index.map(ref[allele1_col])
    if chr_col is not None:
        annotated["CHR"] = variant_index.map(ref[chr_col])
    if bp_col is not None:
        annotated["bp"] = variant_index.map(ref[bp_col])

    if annotated["ID"].isnull().any():
        raise SystemExit("{} results contain unmapped variant_index values".format(label))

    preferred_columns = [
        "variant_index", "ID", "CHR", "bp", "hase_beta_allele", "hase_other_allele"
    ]
    ordered_columns = [col for col in preferred_columns if col in annotated.columns]
    ordered_columns.extend([col for col in annotated.columns if col not in ordered_columns])
    return annotated[ordered_columns]


cohort_pattern = os.path.join(
    run_out, "cohort", "cohort={}".format(study_name), "phenotype={}".format(cpg), "file_*.feather"
)
meta_pattern = os.path.join(
    run_out, "meta", "phenotype={}".format(cpg), "file_*.feather"
)

cohort_df, cohort_files = read_feathers(cohort_pattern, "cohort")
meta_df, meta_files = read_feathers(meta_pattern, "meta")
ref, id_col, allele1_col, allele2_col, chr_col, bp_col = load_reference(reference_file)

cohort_df = annotate_variants(cohort_df, ref, id_col, allele1_col, allele2_col, chr_col, bp_col, "cohort")
meta_df = annotate_variants(meta_df, ref, id_col, allele1_col, allele2_col, chr_col, bp_col, "meta")

cohort_csv = os.path.join(validation_out, "cohort_{}_{}.csv.gz".format(study_name, cpg))
meta_csv = os.path.join(validation_out, "meta_{}.csv.gz".format(cpg))
cohort_df.to_csv(cohort_csv, index=False, compression="gzip")
meta_df.to_csv(meta_csv, index=False, compression="gzip")

print("Combined cohort feather files: {}".format(len(cohort_files)))
print("Combined meta feather files: {}".format(len(meta_files)))
print("Cohort rows: {}".format(cohort_df.shape[0]))
print("Meta rows: {}".format(meta_df.shape[0]))
print("Mapped variant_index using reference: {}".format(reference_file))
print("Wrote cohort CSV: {}".format(cohort_csv))
print("Wrote meta CSV: {}".format(meta_csv))
PY

echo "04e meta-classic validation successfully completed"

echo "Comparing HASE validation output against PLINK positive-control GWAS"

hase_validation_csv="${validation_out}/cohort_${study_name}_${validation_cpg}.csv.gz"
hase_plink_prefix="${study_name}_${validation_cpg}"

check_file "${hase_validation_csv}"

${R_directory}Rscript resources/genetics/plot_hase_vs_plink_validation.R \
    "${hase_validation_csv}" \
    "${plink_positive_control}" \
    "${validation_out}" \
    "${hase_plink_prefix}"

echo "HASE vs PLINK validation plots successfully completed"
