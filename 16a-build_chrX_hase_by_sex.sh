#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_16_dir}/logs_a"
exec &> >(tee ${section_16a_logfile})
print_version

mkdir -p ${hase16_chrx_in_female}
mkdir -p ${hase16_chrx_in_male}
mkdir -p ${hase16_chrx_converting_female}
mkdir -p ${hase16_chrx_converting_male}
mkdir -p ${hase16_chrx_mapping_female}
mkdir -p ${hase16_chrx_mapping_male}
mkdir -p ${light_hase}/data

rm -rf ${hase16_chrx_in_female:?}/*
rm -rf ${hase16_chrx_in_male:?}/*
rm -rf ${hase16_chrx_converting_female:?}/*
rm -rf ${hase16_chrx_converting_male:?}/*
rm -rf ${hase16_chrx_mapping_female:?}/*
rm -rf ${hase16_chrx_mapping_male:?}/*

for ref_file in ref-hrc.ref.gz ref-hrc.ref_info.h5
do
    if [ ! -f "${light_hase}/data/${ref_file}" ]
    then
        if [ ! -f "${hase}/data/${ref_file}" ]
        then
            echo "ERROR: Missing source reference file: ${hase}/data/${ref_file}"
            echo "Please download the HASE reference files from the project SFTP and place them in: ${hase}/data"
            echo "Required files: ref-hrc.ref.gz and ref-hrc.ref_info.h5"
            exit 1
        fi
        echo "Copying reference file: ${ref_file}"
        cp "${hase}/data/${ref_file}" "${light_hase}/data/${ref_file}"
    else
        echo "Reference file already exists: ${light_hase}/data/${ref_file}"
    fi
done

hrc_ref_allele="${light_hase}/data/hrc_ref_allele_16.txt"

zcat ${light_hase}/data/ref-hrc.ref.gz \
    | awk 'NR>1 {print $1 "\t" $4}' \
    > ${hrc_ref_allele}

if [ ! -f "${transformed_methylation_adjusted_pcs}.csv" ]
then
    echo "ERROR: Missing all-probe methylation file: ${transformed_methylation_adjusted_pcs}.csv"
    exit 1
fi

if [ ! -f "${covariates_combined}.txt" ]
then
    echo "ERROR: Missing covariates file: ${covariates_combined}.txt"
    exit 1
fi

if [ ! -f "${bfile}.fam" ]
then
    echo "ERROR: Missing PLINK fam file: ${bfile}.fam"
    exit 1
fi

sex_col=$(awk 'NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "Sex_factor") {
            print i
            exit
        }
    }
}' "${covariates_combined}.txt")

if [ -z "${sex_col}" ]
then
    echo "ERROR: Cannot find Sex_factor column in ${covariates_combined}.txt"
    exit 1
fi

echo "Using Sex_factor column: ${sex_col}"

check_chrx_coding() {
    bim_file="$1"

    if [ ! -f "${bim_file}" ]
    then
        echo "ERROR: Missing BIM file for chrX coding check: ${bim_file}"
        exit 1
    fi

    nX=$(awk '{
        chr=toupper($1)
        id=toupper($2)
        if (chr == "X" || id ~ /^X:/) {
            n++
        }
    } END {print n + 0}' "${bim_file}")

    if [ "$nX" -gt "0" ]
    then
        echo "ERROR: wrong chrX coding in ${bim_file}"
        echo "Found ${nX} rows where the first BIM column is X or the second BIM column starts with X:"
        awk 'BEGIN {OFS="\t"} {
            chr=toupper($1)
            id=toupper($2)
            if (chr == "X" || id ~ /^X:/) {
                print $1, $2
            }
        }' "${bim_file}" | head
        exit 1
    fi
}

make_sex_keep_file() {
    sex_label="$1"
    sex_code="$2"
    id_file="$3"
    keep_file="$4"

    awk -v sex_col="${sex_col}" -v sex_code="${sex_code}" 'NR > 1 && $sex_col == sex_code {print $1}' \
        "${covariates_combined}.txt" > "${id_file}"

    sex_id_count=$(wc -l < "${id_file}" | awk '{print $1}')
    echo "Found ${sex_id_count} ${sex_label} samples in ${covariates_combined}.txt"

    if [ "${sex_id_count}" -eq "0" ]
    then
        echo "Skipping ${sex_label}: no sample IDs found using Sex_factor == ${sex_code}"
        rm -f "${id_file}"
        return 1
    fi

    awk 'BEGIN {OFS="\t"} NR==FNR {ids[$1]; next} $2 in ids {print $1, $2}' \
        "${id_file}" "${bfile}.fam" > "${keep_file}"

    sex_keep_count=$(wc -l < "${keep_file}" | awk '{print $1}')
    echo "Found ${sex_keep_count} ${sex_label} samples in ${bfile}.fam"

    rm -f "${id_file}"

    if [ "${sex_keep_count}" -eq "0" ]
    then
        echo "Skipping ${sex_label}: no sample IDs from ${covariates_combined}.txt were found in ${bfile}.fam"
        rm -f "${keep_file}"
        return 1
    fi

    if [ "${sex_keep_count}" -lt "${sex_id_count}" ]
    then
        echo "WARNING: ${sex_label} samples in fam (${sex_keep_count}) are fewer than IDs in covariates (${sex_id_count})"
    fi

    return 0
}

make_chrx_hase_input() {
    sex_label="$1"
    keep_file="$2"
    sex_input_dir="$3"
    sex_haseinput_pgen="${sex_input_dir}/data_haseinput_pgen"
    sex_bfile_prefix="${sex_input_dir}/data_${sex_label}"

    echo "Preparing chrX genotype data for ${sex_label} samples"

    rm -f "${sex_input_dir}"/*.bed "${sex_input_dir}"/*.bim "${sex_input_dir}"/*.fam "${sex_input_dir}"/*.log "${sex_input_dir}"/*.nosex
    rm -f "${sex_input_dir}"/*.pgen "${sex_input_dir}"/*.pvar "${sex_input_dir}"/*.psam

    ${plink2} \
        --bfile "${bfile}" \
        --keep "${keep_file}" \
        --chr X \
        --sort-vars \
        --set-all-var-ids @:#_\$1_\$2 \
        --ref-allele force "${hrc_ref_allele}" 2 1 \
        --make-pgen \
        --output-chr 26 \
        --out "${sex_haseinput_pgen}" \
        --threads "${nthreads}"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: PLINK2 chrX pgen preparation failed for ${sex_label} samples"
        exit 1
    fi

    ${plink2} \
        --pfile "${sex_haseinput_pgen}" \
        --make-bed \
        --output-chr 26 \
        --out "${sex_bfile_prefix}" \
        --threads "${nthreads}"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: PLINK2 chrX bed conversion failed for ${sex_label} samples"
        exit 1
    fi

    rm -f "${sex_haseinput_pgen}.pgen" "${sex_haseinput_pgen}.pvar" "${sex_haseinput_pgen}.psam" "${sex_haseinput_pgen}.log"

    check_chrx_coding "${sex_bfile_prefix}.bim"

    rm -f "${sex_bfile_prefix}.log"
    rm -f "${keep_file}"
}

convert_and_map_chrx_hase() {
    sex_label="$1"
    sex_input_dir="$2"
    sex_converting_dir="$3"
    sex_mapping_dir="$4"

    echo "Start converting chrX genetic data of ${sex_label} samples"
    python ${light_hase}/hase.py \
        -mode converting \
        -g ${sex_input_dir} \
        -o ${sex_converting_dir} \
        -study_name ${study_name}
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase converting failed for ${sex_label} chrX samples"
        exit 1
    fi

    echo "Start mapping chrX genetic data of ${sex_label} samples"
    python ${light_hase}/tools/mapper.py \
        -g ${sex_converting_dir} \
        -o ${sex_mapping_dir} \
        -study_name ${study_name} \
        -ref_name "ref-hrc"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase mapper failed for ${sex_label} chrX samples"
        exit 1
    fi
}

processed_sex_count=0

if make_sex_keep_file "female" "F" "${hase16_chrx_in_female}/female_id" "${hase16_chrx_in_female}/female_fid_id"
then
    make_chrx_hase_input "female" "${hase16_chrx_in_female}/female_fid_id" "${hase16_chrx_in_female}"
    convert_and_map_chrx_hase "female" "${hase16_chrx_in_female}" "${hase16_chrx_converting_female}" "${hase16_chrx_mapping_female}"
    processed_sex_count=$((processed_sex_count + 1))
fi

if make_sex_keep_file "male" "M" "${hase16_chrx_in_male}/male_id" "${hase16_chrx_in_male}/male_fid_id"
then
    make_chrx_hase_input "male" "${hase16_chrx_in_male}/male_fid_id" "${hase16_chrx_in_male}"
    convert_and_map_chrx_hase "male" "${hase16_chrx_in_male}" "${hase16_chrx_converting_male}" "${hase16_chrx_mapping_male}"
    processed_sex_count=$((processed_sex_count + 1))
fi

if [ "${processed_sex_count}" -eq "0" ]
then
    echo "ERROR: No female or male chrX genotype inputs were prepared for module 16"
    exit 1
fi

echo "Successfully prepared module 16 chrX HASE inputs"
