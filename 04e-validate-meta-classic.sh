#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

validation_cpg="${positive_control_cpg:-cg07959070}"
validation_out="${section_04_dir}/meta_classic_validation_positive_control"
meta_inputs="${section_04_dir}/meta_inputs"
run_out="${validation_out}/run"
selected_covariates="${validation_out}/selected_covariates.tsv"
ph_id_inc="${validation_out}/positive_control_cpg.txt"

mkdir -p "${validation_out}" "${run_out}"
exec &> >(tee "${validation_out}/log.txt")
print_version

echo "Validating 04e meta inputs with light_hase meta-classic"
echo "Study: ${study_name}"
echo "Positive control CpG: ${validation_cpg}"
echo "Meta inputs: ${meta_inputs}"
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

printf "ID\n%s\n" "${validation_cpg}" > "${ph_id_inc}"

if [ -n "${APPTAINER_BIN:-}" ] && [ -n "${HASE_SIF:-}" ]; then
    apptainer_bind="${APPTAINER_BIND:-${home_directory},${scripts_directory}}"
    PYTHON_RUNNER=("${APPTAINER_BIN}" exec)
    if [ -n "${apptainer_bind}" ]; then
        PYTHON_RUNNER+=(--bind "${apptainer_bind}")
    fi
    PYTHON_RUNNER+=("${HASE_SIF}" python)
else
    PYTHON_RUNNER=("${Python_directory}python")
fi

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
    "${validation_cpg}" <<'PY'
import glob
import os
import sys

import pandas as pd

run_out, validation_out, study_name, cpg = sys.argv[1:5]


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


cohort_pattern = os.path.join(
    run_out, "cohort", "cohort={}".format(study_name), "phenotype={}".format(cpg), "file_*.feather"
)
meta_pattern = os.path.join(
    run_out, "meta", "phenotype={}".format(cpg), "file_*.feather"
)

cohort_df, cohort_files = read_feathers(cohort_pattern, "cohort")
meta_df, meta_files = read_feathers(meta_pattern, "meta")

cohort_csv = os.path.join(validation_out, "cohort_{}_{}.csv.gz".format(study_name, cpg))
meta_csv = os.path.join(validation_out, "meta_{}.csv.gz".format(cpg))
cohort_df.to_csv(cohort_csv, index=False, compression="gzip")
meta_df.to_csv(meta_csv, index=False, compression="gzip")

print("Combined cohort feather files: {}".format(len(cohort_files)))
print("Combined meta feather files: {}".format(len(meta_files)))
print("Cohort rows: {}".format(cohort_df.shape[0]))
print("Meta rows: {}".format(meta_df.shape[0]))
print("Wrote cohort CSV: {}".format(cohort_csv))
print("Wrote meta CSV: {}".format(meta_csv))
PY

echo "04e meta-classic validation successfully completed"
