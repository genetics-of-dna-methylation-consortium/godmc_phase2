#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_05a_logfile})
print_version

# this script is to observe the female/male genotype data, and generate mapper

mkdir -p ${hase_in_female}
mkdir -p ${hase_in_male}
mkdir -p ${hase_converting_female}
mkdir -p ${hase_converting_male}
mkdir -p ${hase_mapping_female}
mkdir -p ${hase_mapping_male}
mkdir -p ${light_hase}/data

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

echo "Flipping alleles into hrc reference allele order"

hrc_ref_allele="${light_hase}/data/hrc_ref_allele.txt"

zcat ${light_hase}/data/ref-hrc.ref.gz \
    | awk 'NR>1 {print $1 "\t" $4}' \
    > ${hrc_ref_allele}

echo "ref-hrc.ref.gz lines including header:"
zcat ${light_hase}/data/ref-hrc.ref.gz | wc -l

echo "${hrc_ref_allele} lines without header:"
wc -l ${hrc_ref_allele}

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

check_chr_x_coding() {
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
        echo "First affected rows, showing BIM columns 1 and 2:"
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

# female samples
if [ -f ${transformed_methylation_adjusted_pcs}.Female.chrX.csv ];
then
    awk -v sex_col="${sex_col}" 'NR > 1 && $sex_col == "F" {print $1}' \
        "${covariates_combined}.txt" > "${hase_in_female}/female_id"

    female_id_count=$(wc -l < "${hase_in_female}/female_id" | awk '{print $1}')
    echo "Found ${female_id_count} female samples in ${covariates_combined}.txt"

    if [ "${female_id_count}" -eq "0" ]
    then
        echo "ERROR: No female sample IDs found using Sex_factor == F"
        exit 1
    fi

    awk 'BEGIN {OFS="\t"} NR==FNR {ids[$1]; next} $2 in ids {print $1, $2}' \
        "${hase_in_female}/female_id" "${bfile}.fam" > "${hase_in_female}/female_fid_id"

    female_keep_count=$(wc -l < "${hase_in_female}/female_fid_id" | awk '{print $1}')
    echo "Found ${female_keep_count} female samples in ${bfile}.fam"

    if [ "${female_keep_count}" -eq "0" ]
    then
        echo "ERROR: No female sample IDs from ${covariates_combined}.txt were found in ${bfile}.fam"
        exit 1
    fi

    if [ "${female_keep_count}" -lt "${female_id_count}" ]
    then
        echo "WARNING: female samples in fam (${female_keep_count}) are fewer than IDs in covariates (${female_id_count})"
    fi

    echo "Preparing female genotype data in hrc reference allele order"

    female_haseinput_pgen="${hase_in_female}/data_haseinput_pgen"

    rm -f "${hase_in_female}"/*.bed "${hase_in_female}"/*.bim "${hase_in_female}"/*.fam "${hase_in_female}"/*.log "${hase_in_female}"/*.nosex
    rm -f "${hase_in_female}"/*.pgen "${hase_in_female}"/*.pvar "${hase_in_female}"/*.psam

    ${plink2} \
        --bfile "${bfile}" \
        --keep "${hase_in_female}/female_fid_id" \
        --sort-vars \
        --set-all-var-ids @:#_\$1_\$2 \
        --ref-allele force "${hrc_ref_allele}" 2 1 \
        --make-pgen \
        --output-chr 26 \
        --out "${female_haseinput_pgen}" \
        --threads "${nthreads}"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: PLINK2 pgen preparation failed for female samples"
        exit 1
    fi

    ${plink2} \
        --pfile "${female_haseinput_pgen}" \
        --make-bed \
        --output-chr 26 \
        --out "${hase_in_female}/data" \
        --threads "${nthreads}"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: PLINK2 bed conversion failed for female samples"
        exit 1
    fi

    rm -f "${female_haseinput_pgen}.pgen" "${female_haseinput_pgen}.pvar" "${female_haseinput_pgen}.psam"

    check_chr_x_coding "${hase_in_female}/data.bim"

    rm -f "${hase_in_female}/female_id" "${hase_in_female}/female_fid_id"

    echo "Start converting genetic data of female samples"
    python ${light_hase}/hase.py \
        -mode converting \
        -g ${hase_in_female} \
        -o ${hase_converting_female} \
        -study_name ${study_name} # the name for your study
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase converting failed for female samples"
        exit 1
    fi

    echo "Start mapping genetic data of female samples"
    python ${light_hase}/tools/mapper.py \
        -g ${hase_converting_female} \
        -o ${hase_mapping_female} \
        -study_name ${study_name} \
        -ref_name "ref-hrc"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase mapper failed for female samples"
        exit 1
    fi
else
    echo "file ${transformed_methylation_adjusted_pcs}.Female.chrX.csv does not exist, please check if no female samples in your dataset"
fi

# male samples
if [ -f ${transformed_methylation_adjusted_pcs}.Male.chrX.csv ];
then
    awk -v sex_col="${sex_col}" 'NR > 1 && $sex_col == "M" {print $1}' \
        "${covariates_combined}.txt" > "${hase_in_male}/male_id"

    male_id_count=$(wc -l < "${hase_in_male}/male_id" | awk '{print $1}')
    echo "Found ${male_id_count} male samples in ${covariates_combined}.txt"

    if [ "${male_id_count}" -eq "0" ]
    then
        echo "ERROR: No male sample IDs found using Sex_factor == M"
        exit 1
    fi

    awk 'BEGIN {OFS="\t"} NR==FNR {ids[$1]; next} $2 in ids {print $1, $2}' \
        "${hase_in_male}/male_id" "${bfile}.fam" > "${hase_in_male}/male_fid_id"

    male_keep_count=$(wc -l < "${hase_in_male}/male_fid_id" | awk '{print $1}')
    echo "Found ${male_keep_count} male samples in ${bfile}.fam"

    if [ "${male_keep_count}" -eq "0" ]
    then
        echo "ERROR: No male sample IDs from ${covariates_combined}.txt were found in ${bfile}.fam"
        exit 1
    fi

    if [ "${male_keep_count}" -lt "${male_id_count}" ]
    then
        echo "WARNING: male samples in fam (${male_keep_count}) are fewer than IDs in covariates (${male_id_count})"
    fi

    echo "Preparing male genotype data in hrc reference allele order"

    male_haseinput_pgen="${hase_in_male}/data_haseinput_pgen"

    rm -f "${hase_in_male}"/*.bed "${hase_in_male}"/*.bim "${hase_in_male}"/*.fam "${hase_in_male}"/*.log "${hase_in_male}"/*.nosex
    rm -f "${hase_in_male}"/*.pgen "${hase_in_male}"/*.pvar "${hase_in_male}"/*.psam

    ${plink2} \
        --bfile "${bfile}" \
        --keep "${hase_in_male}/male_fid_id" \
        --sort-vars \
        --set-all-var-ids @:#_\$1_\$2 \
        --ref-allele force "${hrc_ref_allele}" 2 1 \
        --make-pgen \
        --output-chr 26 \
        --out "${male_haseinput_pgen}" \
        --threads "${nthreads}"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: PLINK2 pgen preparation failed for male samples"
        exit 1
    fi

    ${plink2} \
        --pfile "${male_haseinput_pgen}" \
        --make-bed \
        --output-chr 26 \
        --out "${hase_in_male}/data" \
        --threads "${nthreads}"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: PLINK2 bed conversion failed for male samples"
        exit 1
    fi

    rm -f "${male_haseinput_pgen}.pgen" "${male_haseinput_pgen}.pvar" "${male_haseinput_pgen}.psam"

    check_chr_x_coding "${hase_in_male}/data.bim"

    rm -f "${hase_in_male}/male_id" "${hase_in_male}/male_fid_id"

    echo "Start converting genetic data of male samples"
    python ${light_hase}/hase.py \
        -mode converting \
        -g ${hase_in_male} \
        -o ${hase_converting_male} \
        -study_name ${study_name} # the name for your study
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase converting failed for male samples"
        exit 1
    fi

    echo "Start mapping genetic data of male samples"
    python ${light_hase}/tools/mapper.py \
        -g ${hase_converting_male} \
        -o ${hase_mapping_male} \
        -study_name ${study_name} \
        -ref_name "ref-hrc"
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase mapper failed for male samples"
        exit 1
    fi
else
    echo "file ${transformed_methylation_adjusted_pcs}.Male.chrX.csv does not exist, please check if no male samples in your dataset"
fi

echo "Successfully finished 05a script"
