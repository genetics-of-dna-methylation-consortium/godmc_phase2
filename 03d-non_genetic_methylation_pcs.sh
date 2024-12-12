#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_03d_logfile})
print_version

if [ -f ${meth_pcs_transformed}.txt ]; then
echo "Performing genetic analysis of methylation PCs of transformed data by MatrixEQTL"
${R_directory}Rscript resources/methylation/genetic_meth_pcs.R \
	${tabfile}.tab \
	${meth_pcs_transformed} \
	${nongenetic_meth_pcs_transformed} \
	${meth_pc_threshold} \
	${nthreads}
fi
    

if [ -f ${meth_pcs_untransformed}.txt ]; then
echo "Performing genetic analysis of methylation PCs of untransformed data by MatrixEQTL"
${R_directory}Rscript resources/methylation/genetic_meth_pcs.R \
	${tabfile}.tab \
	${meth_pcs_untransformed} \
	${nongenetic_meth_pcs_untransformed} \
	${meth_pc_threshold} \
	${nthreads}
fi
echo "Successfully generating non-genetic methylation PCs"
