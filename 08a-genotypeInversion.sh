#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_08a_logfile})
print_version

# Change annotation ####
#'#################################################################################

echo "Adapting plink files"

mkdir -p ${inv_processed_dir}
## Remove indels
${plink2} --bfile  ${bfile} --extract range resources/inversions/inversion_ranges.txt --make-bed --out ${inv_processed_dir}/invSNPs

echo "Inferring inversions"

${R_directory}Rscript --vanilla resources/inversions/genotypeInversions.R ${inv_processed_dir}/invSNPs ${inv_processed_dir} ${section_08_dir} ${nthreads} ${inversion_chunks}

echo "Successfully inferred inversions"


