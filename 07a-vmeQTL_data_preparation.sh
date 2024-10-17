#!/bin/bash

./resources/setup.sh
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
