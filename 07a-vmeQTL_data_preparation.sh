#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_07a_logfile})
print_version

# create a directory to store vmeQTL processed result
echo "check the existence of ${meth_vmeQTL_directory}"

if [ ! -d ${meth_vmeQTL_directory} ]
then
    echo "${meth_vmeQTL_directory} does not exist, make this directory"
    mkdir ${meth_vmeQTL_directory}
fi

echo "Split methylation data into each chromosome"

if [ "$genetic_chunks" -gt 100 ]
then
    echo "#####################################################################################################################"
    echo "#####################################################################################################################"
    echo "#####################################################################################################################"
    echo "WARNING: The genetic chunk in config file is ${genetic_chunk}"
    echo "We would highly recommend to change genetic_chunk to 100 to save your running time of 07b, unless your sample size is large"
    echo "If you change genetic_chunk for module 07, please rerun 02b"
    echo "Please check the wiki https://github.com/genetics-of-dna-methylation-consortium/godmc_phase2/wiki/Run-variance-meQTL-analysis or contact xiaopu.1.zhang@kcl.ac.uk if you are unsure before you running 07"
    echo "#####################################################################################################################"
    echo "#####################################################################################################################"
    echo "#####################################################################################################################"

fi

${R_directory}Rscript \
    ./resources/methylation/vmeQTL_process_tabfile.R \
    ${untransformed_methylation_adjusted_pcs}.RData \
    ${meth_vmeQTL_input_chr}

for chr in $(seq 1 22)
do
    echo "convert chr ${chr} methylation data to bod format"
    ${osca} \
        --tefile ${meth_vmeQTL_input_chr}${chr} \
        --methylation-m \
        --make-bod \
        --no-fid \
        --out ${meth_vmeQTL_input_chr}${chr}

    echo "update annotation of chr ${chr} methylation data"
    ${osca} \
        --befile ${meth_vmeQTL_input_chr}${chr} \
        --update-opi ${meth_vmeQTL_annotation}.opi

done

echo "Convert chunked genetic data to plink data format"

for chunk in $(seq 1 ${genetic_chunks})
do

startsnp=`cut -f 1 ${tabfile}.tab.${chunk} | grep -v snpid | head -n1`
endsnp=`cut -f 1 ${tabfile}.tab.${chunk} | tail -n1`

${plink} \
--bfile ${bfile} \
--snps ${startsnp}-${endsnp} \
--make-bed \
--out ${tabfile}.tab.${chunk}

done

echo "Data preparation for vmeQTL analysis is successful"
