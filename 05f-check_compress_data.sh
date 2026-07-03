#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_05_dir}/logs_f"
exec &> >(tee "${section_05f_logfile}")
print_version

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

check_any_file() {
    dir="$1"
    pattern="$2"
    label="$3"

    if ! find "${dir}" -maxdepth 1 -type f -name "${pattern}" | grep -q .; then
        fail "Missing ${label} in ${dir}"
    fi
}

check_meta_inputs() {
    sex_label="$1"
    meta_inputs="${section_05_dir}/meta_inputs_${sex_label}"

    echo "Checking Module 05 ${sex_label} meta-analysis inputs: ${meta_inputs}"

    check_dir "${meta_inputs}"
    check_dir "${meta_inputs}/part_dev"
    check_dir "${meta_inputs}/mapping"
    check_dir "${meta_inputs}/use_data"
    check_dir "${meta_inputs}/use_data/genotype"
    check_dir "${meta_inputs}/use_data/individuals"
    check_dir "${meta_inputs}/use_data/probes"
    check_dir "${meta_inputs}/use_data/phenotypes"

    check_any_file "${meta_inputs}/part_dev" "*.npy" "single-site partial derivative npy files"
    check_any_file "${meta_inputs}/mapping" "*.npy" "mapper npy files"
    check_any_file "${meta_inputs}/use_data/genotype" "*.h5" "encoded genotype h5 files"
    check_any_file "${meta_inputs}/use_data/individuals" "*.h5" "encoded individual h5 files"
    check_any_file "${meta_inputs}/use_data/probes" "*.h5" "probe h5 files"
    check_any_file "${meta_inputs}/use_data/phenotypes" "*.csv" "encoded phenotype csv files"
}

check_positive_control_outputs() {
    sex_label="$1"
    validation_cpg="${sexchr_positive_control_cpg}"

    echo "Checking Module 05 ${sex_label} positive-control validation outputs for ${validation_cpg}"

    check_file "${section_05_dir}/sexchr_positive_control_validation/hase/${sex_label}/cohort_${study_name}_${validation_cpg}.csv.gz"
    check_file "${section_05_dir}/sexchr_positive_control_validation/hase/${sex_label}/meta_${validation_cpg}.csv.gz"
    check_file "${section_05_dir}/sexchr_positive_control_validation/plink/${sex_label}/positive_control_${sex_label}_${validation_cpg}.PHENO1.glm.linear.gz"
    check_file "${section_05_dir}/sexchr_positive_control_validation/hase/${sex_label}/${study_name}_${sex_label}_${validation_cpg}.merged.tsv.gz"
}

check_expected_sex_outputs() {
    sex_label="$1"

    check_meta_inputs "${sex_label}"
    check_positive_control_outputs "${sex_label}"
}

if [ "${config_file:0:1}" = "/" ]; then
    config_to_archive="${config_file}"
else
    config_to_archive="${scripts_directory}/${config_file}"
fi

check_file "${config_to_archive}"
check_file "${scripts_directory}/resources/parameters"
check_file "${covariates_combined}.txt"

sex_col=$(awk 'NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "Sex_factor") {
            print i
            exit
        }
    }
}' "${covariates_combined}.txt")

if [ -z "${sex_col}" ]; then
    fail "Cannot find Sex_factor column in ${covariates_combined}.txt"
fi

n_female=$(awk -v sex_col="${sex_col}" 'NR > 1 && $sex_col == "F" {n++} END {print n + 0}' "${covariates_combined}.txt")
n_male=$(awk -v sex_col="${sex_col}" 'NR > 1 && $sex_col == "M" {n++} END {print n + 0}' "${covariates_combined}.txt")

echo "Sex_factor counts: female=${n_female}, male=${n_male}"

if [ "${n_female}" -gt "0" ] && [ "${n_male}" -gt "0" ]; then
    echo "Cohort contains both female and male samples"
    check_expected_sex_outputs "female"
    check_expected_sex_outputs "male"
elif [ "${n_female}" -gt "0" ]; then
    echo "Cohort female only"
    check_expected_sex_outputs "female"
elif [ "${n_male}" -gt "0" ]; then
    echo "Cohort male only"
    check_expected_sex_outputs "male"
else
    fail "No M or F values found in Sex_factor column of ${covariates_combined}.txt"
fi

cd "${home_directory}"

archive="${home_directory}/results/05_${study_name}.tgz"
checksum="${archive}.md5sum"
encrypted="${archive}.gpg"

echo "Removing old tgz, md5 and gpg files"
rm -f "${archive}" "${checksum}" "${encrypted}"

echo "Compressing Module 05 results"
tar -zcf "${archive}" \
    --exclude="results/05/sexchr_positive_control_validation/hase/*/run/*"  \
    "${config_to_archive}" \
    "${scripts_directory}/resources/parameters" \
    "results/05"

check_file "${archive}"

echo "Generating md5 checksum"
cd "${home_directory}/results"
md5sum "05_${study_name}.tgz" > "05_${study_name}.tgz.md5sum"
md5sum -c "05_${study_name}.tgz.md5sum"

# echo "Encrypting Module 05 archive"
# gpg --output "05_${study_name}.tgz.gpg" \
#     --symmetric \
#     --cipher-algo AES256 \
#     "05_${study_name}.tgz"

# check_file "${checksum}"
# check_file "${encrypted}"

echo ""
echo "Module 05 archive successfully created and encrypted."
echo "Please upload these files to Google Drive: https://drive.google.com/drive/folders/1q1djBVG5ms-Ud3btZmIf_YHHnjgphxem?usp=share_link"
echo "1. ${checksum}"
echo "2. ${encrypted}"
