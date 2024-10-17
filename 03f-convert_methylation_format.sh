#!/bin/bash

./resources/setup.sh "$@"
exec &> >(tee ${section_03f_logfile})
print_version

mkdir -p ${hase_cov}
mkdir -p ${hase_cov_males}
mkdir -p ${hase_cov_females}

${R_directory}Rscript resources/methylation/methylation_hase_format.R \
		${untransformed_methylation_adjusted_pcs} \
		${transformed_methylation_adjusted_pcs} \
		${covariates_combined}.txt \
		${hase_cov} \
                ${hase_cov_males} \
                ${hase_cov_females}

echo "Successfully converted methylation data"

