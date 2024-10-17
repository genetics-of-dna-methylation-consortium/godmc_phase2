#!/bin/bash

./resources/setup.sh "$@"
mkdir -p ${section_08_dir}
mkdir -p ${section_08_dir}/logs_b/
set -e
exec &> >(tee ${section_08b_logfile})
print_version

geno="${inv_processed_dir}/inversionQTL.txt"
methy="${transformed_methylation_adjusted_pcs}.RData"
methy_trans="${inv_processed_dir}/methy_transposed.txt"

out="${section_08_dir}/invmeqtl.Rdata"
threshold=1

echo "Transposing methylation matrix for MatrixEQTL"
Rscript resources/inversions/adapt_methylation_meQTL.R ${methy} ${methy_trans}

echo "Performing inversion meQTL analysis"

${R_directory}Rscript resources/genetics/matrixeqtl.R ${geno} ${methy_trans} NULL ${threshold} ${out}


echo "Successfully completed inversion meQTL analysis"
