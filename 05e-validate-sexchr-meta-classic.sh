#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_05_dir}/logs_e"
exec &> >(tee "${section_05e_logfile}")
print_version

validation_cpg="${sexchr_positive_control_cpg}"
validation_out="${section_05_dir}/sexchr_positive_control_validation"
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

if [ -z "${validation_cpg}" ]; then
    fail "sexchr_positive_control_cpg is empty. Please set it in your config before running 05e."
fi

mkdir -p "${plink_out}/female" "${plink_out}/male"
sex_final_validation_count=0

echo "Running 05e final sex-stratified PLINK/R validation"
echo "Study: ${study_name}"
echo "Positive control CpG: ${validation_cpg}"
echo "Output: ${validation_out}"

run_sex_final_validation() {
    sex_label="$1"
    pheno_dir="$2"
    bfile_prefix="$3"
    out_dir="$4"
    hase_dir="$5"

    phenotype_csv="${pheno_dir}/methylation_data.csv"
    extracted_csv="${out_dir}/${validation_cpg}.positive_control.csv"
    plink_pheno="${out_dir}/${validation_cpg}.positive_control.plink"
    plink_prefix="${out_dir}/positive_control_${sex_label}_${validation_cpg}"
    plink_glm="${plink_prefix}.PHENO1.glm.linear"
    plink_glm_gz="${plink_glm}.gz"
    plot_file_list="${out_dir}/positive.control.${sex_label}.file.txt"
    hase_validation_csv="${hase_dir}/cohort_${study_name}_${validation_cpg}.csv.gz"
    hase_plink_prefix="${study_name}_${sex_label}_${validation_cpg}"

    echo "Preparing ${sex_label} final PLINK/R validation"
    if [ ! -f "${phenotype_csv}" ]; then
        echo "Skipping ${sex_label} final validation because phenotype input is missing: ${phenotype_csv}"
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

    ${R_directory}Rscript resources/genetics/make_control.R \
        "${extracted_csv}" \
        "${bfile_prefix}.fam" \
        "${plink_pheno}"

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
        "${sexchr_positive_control_snp_chr}" \
        "${sexchr_positive_control_snp_pos}" \
        "${sexchr_positive_control_snp_window}" \
        "${sexchr_positive_control_threshold}"
    if [ "$?" -ne "0" ]; then
        fail "plot_gwas.R failed for ${sex_label} PLINK validation"
    fi

    echo "Comparing ${sex_label} HASE validation output against PLINK positive-control GWAS"
    ${R_directory}Rscript resources/genetics/plot_hase_vs_plink_validation.R \
        "${hase_validation_csv}" \
        "${plink_glm_gz}" \
        "${hase_dir}" \
        "${hase_plink_prefix}"
    if [ "$?" -ne "0" ]; then
        fail "plot_hase_vs_plink_validation.R failed for ${sex_label}"
    fi

    sex_final_validation_count=$((sex_final_validation_count + 1))
}

run_sex_final_validation \
    "female" \
    "${hase_pheno_female}" \
    "${hase_in_female}/data_female" \
    "${plink_out}/female" \
    "${hase_out}/female"

run_sex_final_validation \
    "male" \
    "${hase_pheno_male}" \
    "${hase_in_male}/data_male" \
    "${plink_out}/male" \
    "${hase_out}/male"

if [ "${sex_final_validation_count}" -eq "0" ]; then
    fail "No sex-specific phenotype inputs were found for 05e; nothing to validate."
fi

echo "05e final sex-stratified PLINK/R validation successfully completed"
