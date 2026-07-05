#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_16_dir}/logs_b"
exec &> >(tee ${section_16b_logfile})
print_version

mkdir -p ${hase16_chrx_encoding_female}
mkdir -p ${hase16_chrx_encoding_male}
mkdir -p ${hase16_allprobes_pheno_female}
mkdir -p ${hase16_allprobes_pheno_male}

rm -rf ${hase16_chrx_encoding_female:?}/*
rm -rf ${hase16_chrx_encoding_male:?}/*
rm -rf ${hase16_allprobes_pheno_female:?}/*
rm -rf ${hase16_allprobes_pheno_male:?}/*

if [ ! -f "${transformed_methylation_adjusted_pcs}.csv" ]
then
    echo "ERROR: Missing autosomal methylation file: ${transformed_methylation_adjusted_pcs}.csv"
    exit 1
fi

check_file() {
    if [ ! -f "$1" ]
    then
        echo "ERROR: Missing required file: $1"
        exit 1
    fi
}

append_methylation_csv_by_fam() {
    sex_label="$1"
    fam_file="$2"
    methylation_csv="$3"
    output_csv="$4"
    write_header="$5"

    check_file "${fam_file}"
    check_file "${methylation_csv}"

    fam_count=$(wc -l < "${fam_file}" | awk '{print $1}')
    if [ "${fam_count}" -eq "0" ]
    then
        echo "ERROR: ${sex_label} chrX genotype fam file contains no samples: ${fam_file}"
        exit 1
    fi

    awk -F',' -v OFS=',' -v sex_label="${sex_label}" -v write_header="${write_header}" '
        NR == FNR {
            split($0, fam_fields, /[[:space:]]+/)
            if (fam_fields[2] == "") {
                print "ERROR: Missing IID in fam file line " FNR ": " FILENAME > "/dev/stderr"
                exit 2
            }
            n++
            ids[n] = fam_fields[2]
            next
        }
        FNR == 1 {
            for (i = 1; i <= NF; i++) {
                col[$i] = i
            }
            missing = 0
            for (i = 1; i <= n; i++) {
                if (!(ids[i] in col)) {
                    print "ERROR: Missing " sex_label " sample " ids[i] " in methylation CSV header: " FILENAME > "/dev/stderr"
                    missing = 1
                }
            }
            if (missing) {
                exit 2
            }
            if (write_header == "1") {
                printf "%s", $1
                for (i = 1; i <= n; i++) {
                    printf "%s%s", OFS, ids[i]
                }
                printf "\n"
            }
            next
        }
        {
            printf "%s", $1
            for (i = 1; i <= n; i++) {
                printf "%s%s", OFS, $(col[ids[i]])
            }
            printf "\n"
        }
    ' "${fam_file}" "${methylation_csv}" >> "${output_csv}"

    if [ "$?" -ne "0" ]
    then
        echo "ERROR: Failed to append ${methylation_csv} to ${output_csv} for ${sex_label}"
        exit 1
    fi
}

build_sex_specific_phenotype() {
    sex_label="$1"
    fam_file="$2"
    pheno_dir="$3"
    shift 3

    output_csv="${pheno_dir}/methylation_data.csv"
    rm -f "${output_csv}"

    echo "Building ${sex_label} sex-specific all-probe methylation phenotype: ${output_csv}"

    write_header=1
    for methylation_csv in "$@"
    do
        append_methylation_csv_by_fam \
            "${sex_label}" \
            "${fam_file}" \
            "${methylation_csv}" \
            "${output_csv}" \
            "${write_header}"
        write_header=0
    done
}

encode_chrx_all_probes() {
    sex_label="$1"
    converting_dir="$2"
    mapping_dir="$3"
    fam_file="$4"
    pheno_dir="$5"
    encoding_dir="$6"
    shift 6

    if [ ! -d "${converting_dir}/probes" ]
    then
        echo "Skipping ${sex_label}: chrX converted genotype directory is missing: ${converting_dir}/probes"
        return 1
    fi

    if [ ! -d "${mapping_dir}" ]
    then
        echo "Skipping ${sex_label}: chrX mapping directory is missing: ${mapping_dir}"
        return 1
    fi

    build_sex_specific_phenotype "${sex_label}" "${fam_file}" "${pheno_dir}" "$@"

    echo "Encoding sex-specific all available methylation probes against chrX genotypes for ${sex_label} samples"
    python ${light_hase}/hase.py \
        -mode encoding \
        -study_name ${study_name} \
        -g ${converting_dir} \
        -o ${encoding_dir} \
        -mapper ${mapping_dir} \
        -ph ${pheno_dir} \
        -ref_name ref-hrc
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase encoding failed for ${sex_label} module 16 inputs"
        exit 1
    fi

    return 0
}

encoded_sex_count=0

if encode_chrx_all_probes \
    "female" \
    "${hase16_chrx_converting_female}" \
    "${hase16_chrx_mapping_female}" \
    "${hase16_chrx_in_female}/data_female.fam" \
    "${hase16_allprobes_pheno_female}" \
    "${hase16_chrx_encoding_female}" \
    "${transformed_methylation_adjusted_pcs}.csv" \
    "${transformed_methylation_adjusted_pcs}.Female.chrX.csv"
then
    encoded_sex_count=$((encoded_sex_count + 1))
fi

if encode_chrx_all_probes \
    "male" \
    "${hase16_chrx_converting_male}" \
    "${hase16_chrx_mapping_male}" \
    "${hase16_chrx_in_male}/data_male.fam" \
    "${hase16_allprobes_pheno_male}" \
    "${hase16_chrx_encoding_male}" \
    "${transformed_methylation_adjusted_pcs}.csv" \
    "${transformed_methylation_adjusted_pcs}.Male.chrX.csv" \
    "${transformed_methylation_adjusted_pcs}.Male.chrY.csv"
then
    encoded_sex_count=$((encoded_sex_count + 1))
fi

if [ "${encoded_sex_count}" -eq "0" ]
then
    echo "ERROR: No female or male module 16 inputs were encoded"
    exit 1
fi

echo "Successfully encoded module 16 chrX genotype and sex-specific all-probe methylation data"
