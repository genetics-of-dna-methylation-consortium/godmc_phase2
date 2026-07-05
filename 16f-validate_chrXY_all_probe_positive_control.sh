#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_16_dir}/logs_f"
exec &> >(tee "${section_16f_logfile}")
print_version

validation_cpg="${module16_positive_control_cpg}"
validation_out="${section_16_dir}/positive_control_validation"
hase_out="${validation_out}/hase"
plink_out="${validation_out}/plink"

fail() {
    echo "ERROR: $*" >&2
    exit 1
}

check_file() {
    if [ ! -f "$1" ]; then
        fail "Missing required file: $1"
    fi
}

check_bfile() {
    prefix="$1"
    check_file "${prefix}.bed"
    check_file "${prefix}.bim"
    check_file "${prefix}.fam"
}

cpg_exists_in_phenotype() {
    phenotype_csv="$1"
    cpg="$2"

    awk -F',' -v cpg="${cpg}" 'NR > 1 && $1 == cpg {found = 1; exit} END {exit found ? 0 : 1}' \
        "${phenotype_csv}"
}

if [ -z "${validation_cpg}" ]; then
    fail "module16_positive_control_cpg is empty. Please set it in your config before running 16f."
fi

mkdir -p "${plink_out}/female" "${plink_out}/male"
sex_final_validation_count=0

echo "Running 16f Module 16 PLINK/R positive-control validation"
echo "Study: ${study_name}"
echo "Positive control CpG: ${validation_cpg}"
echo "Output: ${validation_out}"
echo "This script uses the R/PLINK environment."

run_sex_final_validation() {
    sex_label="$1"
    pheno_dir="$2"
    bfile_prefix="$3"
    out_dir="$4"
    hase_dir="$5"

    phenotype_csv="${pheno_dir}/methylation_data.csv"
    positive_control_prefix="${transformed_methylation_adjusted}.module16_${sex_label}_${validation_cpg}.positive_control"
    extracted_csv="${positive_control_prefix}"
    plink_ids="${positive_control_prefix}.ids"
    plink_pheno="${positive_control_prefix}.plink"
    plink_prefix="${out_dir}/positive_control_${sex_label}_${validation_cpg}"
    plink_glm="${plink_prefix}.PHENO1.glm.linear"
    plink_glm_gz="${plink_glm}.gz"
    plot_file_list="${out_dir}/positive.control.${sex_label}.file.txt"
    hase_validation_csv="${hase_dir}/cohort_${study_name}_${validation_cpg}.csv.gz"
    hase_plink_prefix="${study_name}_${sex_label}_${validation_cpg}"

    echo "Preparing ${sex_label} Module 16 PLINK/R positive-control validation"
    if [ ! -f "${phenotype_csv}" ]; then
        echo "Skipping ${sex_label}: phenotype input is missing: ${phenotype_csv}"
        return 0
    fi
    if ! cpg_exists_in_phenotype "${phenotype_csv}" "${validation_cpg}"; then
        echo "Skipping ${sex_label}: positive control CpG ${validation_cpg} was not found in ${phenotype_csv}"
        return 0
    fi

    check_bfile "${bfile_prefix}"
    check_file "${hase_validation_csv}"

    awk -F',' -v cpg="${validation_cpg}" 'NR == 1 || $1 == cpg {print $0}' \
        "${phenotype_csv}" > "${extracted_csv}"

    nrow=$(wc -l < "${extracted_csv}" | awk '{print $1}')
    if [ "${nrow}" -lt "2" ]; then
        fail "Positive control CpG ${validation_cpg} was not found in ${phenotype_csv}"
    fi

    awk 'BEGIN {OFS="\t"} NR==FNR {ids[$2]; next} $2 in ids {print $1, $2}' \
        "${intersect_ids_plink}" \
        "${bfile_prefix}.fam" > "${plink_ids}"

    nids=$(wc -l < "${plink_ids}" | awk '{print $1}')
    if [ "${nids}" -eq "0" ]; then
        fail "No ${sex_label} samples from ${bfile_prefix}.fam were found in ${intersect_ids_plink}"
    fi

    ${R_directory}Rscript resources/genetics/make_control.R \
        "${extracted_csv}" \
        "${plink_ids}" \
        "${plink_pheno}"
    if [ "$?" -ne "0" ]; then
        fail "make_control.R failed for ${sex_label} Module 16 positive control"
    fi

    echo "Running PLINK2 for ${sex_label}"
    ${plink2} \
        --bfile "${bfile_prefix}" \
        --pheno "${plink_pheno}" \
        --glm allow-no-covars \
        --allow-extra-chr \
        --human \
        --output-chr 26 \
        --threads "${nthreads}" \
        --out "${plink_prefix}"
    if [ "$?" -ne "0" ]; then
        fail "PLINK2 failed for ${sex_label} Module 16 positive control"
    fi

    check_file "${plink_glm}"

    tr -s " " < "${plink_glm}" | gzip -c > "${plink_glm_gz}"
    rm "${plink_glm}"

    check_file "${plink_glm_gz}"
    echo "Wrote ${sex_label} PLINK result: ${plink_glm_gz}"

    echo "Making ${sex_label} Manhattan and QQ plots"
    echo "${plink_glm_gz}" > "${plot_file_list}"
    ${R_directory}Rscript resources/genetics/plot_gwas.R \
        "${plot_file_list}" \
        12 \
        1 \
        2 \
        3 \
        TRUE \
        "${module16_positive_control_snp_chr}" \
        "${module16_positive_control_snp_pos}" \
        "${module16_positive_control_snp_window}" \
        "${module16_positive_control_threshold}"
    if [ "$?" -ne "0" ]; then
        fail "plot_gwas.R failed for ${sex_label} Module 16 positive control"
    fi

    echo "Comparing ${sex_label} HASE validation output against PLINK positive-control GWAS"
    ${R_directory}Rscript resources/genetics/plot_hase_vs_plink_validation.R \
        "${hase_validation_csv}" \
        "${plink_glm_gz}" \
        "${hase_dir}" \
        "${hase_plink_prefix}"
    if [ "$?" -ne "0" ]; then
        fail "plot_hase_vs_plink_validation.R failed for ${sex_label} Module 16 positive control"
    fi

    sex_final_validation_count=$((sex_final_validation_count + 1))
}

run_sex_final_validation \
    "female" \
    "${hase16_allprobes_pheno_female}" \
    "${hase16_chrx_in_female}/data_female" \
    "${plink_out}/female" \
    "${hase_out}/female"

run_sex_final_validation \
    "male" \
    "${hase16_allprobes_pheno_male}" \
    "${hase16_chrx_in_male}/data_male" \
    "${plink_out}/male" \
    "${hase_out}/male"

if [ "${sex_final_validation_count}" -eq "0" ]; then
    fail "No sex-specific Module 16 PLINK/R positive-control validations were completed."
fi

echo "16f Module 16 PLINK/R positive-control validation successfully completed"
