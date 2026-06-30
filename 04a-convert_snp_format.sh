#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04a_logfile})
print_version
#Please read resources/bin/hase/README_2.md
#An expample is also provided below

mkdir -p ${light_hase_dir_in}
mkdir -p ${light_hase_converting}
mkdir -p ${light_hase}/data

# assuming ref is in ori hase folder
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

# Allele-direction note:
# The HASE/light_hase PLINK reader decodes .bed genotypes as .bim allele2
# dosage, so downstream HASE beta is relative to .bim allele2. The PLINK2
# --ref-allele force step below sets the requested HRC allele as PLINK REF/A2
# (.bim column 6), which PLINK association output often treats as the other
# allele rather than the tested A1 allele. With a straight mapper match in 04c,
# HASE beta is relative to ref str_allele2, and ref str_allele1 is the HASE
# other allele.

check_chr_x_coding() {
    bim_file="$1"
    nX=$(awk '{
        chr=toupper($1)
        id=toupper($2)
        if (chr == "X" || id ~ /^X([:_]|$)/) {
            n++
        }
    } END {print n + 0}' "${bim_file}")

    if [ "$nX" -gt "0" ]
    then
        echo "ERROR: wrong chrX coding in ${bim_file}"
        echo "Found ${nX} rows where the first BIM column is X or the second BIM column starts with X: or X_"
        echo "First affected rows, showing BIM columns 1 and 2:"
        awk 'BEGIN {OFS="\t"} {
            chr=toupper($1)
            id=toupper($2)
            if (chr == "X" || id ~ /^X([:_]|$)/) {
                print $1, $2
            }
        }' "${bim_file}" | head
        exit 1
    fi
}

haseinput_pgen="${bfile}_haseinput_pgen"

${plink2} \
    --bfile "${bfile}" \
    --sort-vars \
    --set-all-var-ids @:#_\$1_\$2 \
    --ref-allele force ${hrc_ref_allele} 2 1 \
    --make-pgen \
    --output-chr 26 \
    --out "${haseinput_pgen}" \
	--threads "${nthreads}"

${plink2} \
    --pfile "${haseinput_pgen}" \
    --make-bed \
    --output-chr 26 \
    --out "${bfile}_haseinput" \
	--threads "${nthreads}"

rm -f "${haseinput_pgen}.pgen" "${haseinput_pgen}.pvar" "${haseinput_pgen}.psam" "${haseinput_pgen}.log"

echo "Cleaning up the input files"

rm -f ${light_hase_dir_in}/*.bed ${light_hase_dir_in}/*.bim ${light_hase_dir_in}/*.fam ${light_hase_dir_in}/*.log ${light_hase_dir_in}/*.nosex ${light_hase_dir_in}/*.pgen ${light_hase_dir_in}/*.pvar ${light_hase_dir_in}/*.psam

cp ${bfile}_haseinput.bim ${light_hase_dir_in}/data.bim
cp ${bfile}_haseinput.fam ${light_hase_dir_in}/data.fam
cp ${bfile}_haseinput.bed ${light_hase_dir_in}/data.bed

check_chr_x_coding "${light_hase_dir_in}/data.bim"

python ${light_hase}/hase.py \
    -mode converting \
    -g ${light_hase_dir_in} \
    -o ${light_hase_converting} \
    -study_name ${study_name} # the name for your study 

echo "Successfully converted genetic data"
