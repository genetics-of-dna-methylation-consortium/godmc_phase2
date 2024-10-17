#!/bin/bash

./resources/setup.sh
exec &> >(tee ${section_03c_logfile})
print_version

phenfile="NULL"

echo "Calculating methylation PCs of transformed DNA methylation data"
${R_directory}Rscript resources/methylation/methylation_pcs.R \
    ${transformed_methylation_adjusted}.RData \
    ${meth_pc_cutoff} \
    ${n_meth_pcs} \
    ${phenfile} \
    ${meth_pcs_transformed}

echo "Calculating methylation PCs of untransformed DNA methylation data"
${R_directory}Rscript resources/methylation/methylation_pcs.R \
    ${untransformed_methylation_adjusted}.RData \
    ${meth_pc_cutoff} \
    ${n_meth_pcs} \
    ${phenfile} \
    ${meth_pcs_untransformed}

echo "Successfully generated non-genetic methylation PCs"
