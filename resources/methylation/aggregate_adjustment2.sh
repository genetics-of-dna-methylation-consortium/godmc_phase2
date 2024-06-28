#!/bin/bash

set -e
source config
exec &> >(tee ${section_03e_logfile}_aggregation)

echo "Aggregating methylation_adjustment2 chunks"
${R_directory}Rscript ${scripts_directory}/resources/methylation/aggregate_chunks.R \
	${transformed_methylation_adjusted_pcs} \
	${meth_chunks} \
	${methylation_array} \
    	${methylation_no_outliers}

${R_directory}Rscript ${scripts_directory}/resources/methylation/aggregate_chunks.R \
	${untransformed_methylation_adjusted_pcs} \
	${meth_chunks} \
    	${methylation_array} \
    	${methylation_no_outliers}

if [ -f ${transformed_methylation_adjusted_pcs}.RData ];then
    rm ${transformed_methylation_adjusted_pcs}.[0-9]*.RData
    rm ${transformed_methylation_adjusted}.[0-9]*.RData
fi

if [ -f ${transformed_methylation_adjusted_pcs}.Female.chrX.RData ]; then
    rm ${transformed_methylation_adjusted}.Female.chrX.[0-9]*.RData
fi

if [ -f ${transformed_methylation_adjusted_pcs}.Male.chrX.RData ]; then
    rm ${transformed_methylation_adjusted}.Male.chrX.[0-9]*.RData
fi

if [ -f ${transformed_methylation_adjusted_pcs}.Male.chrY.RData ]; then
    rm ${transformed_methylation_adjusted}.Male.chrY.[0-9]*.RData
fi

if [ -f ${untransformed_methylation_adjusted_pcs}.RData ]; then
    rm ${untransformed_methylation_adjusted_pcs}.[0-9]*.RData
    rm ${untransformed_methylation_adjusted}.[0-9]*.RData
fi

if [ -f ${untransformed_methylation_adjusted_pcs}.Female.chrX.RData ]; then
    rm ${untransformed_methylation_adjusted}.Female.chrX.[0-9]*.RData
fi

if [ -f ${untransformed_methylation_adjusted_pcs}.Male.chrX.RData ]; then
    rm ${untransformed_methylation_adjusted}.Male.chrX.[0-9]*.RData
fi

if [ -f ${untransformed_methylation_adjusted_pcs}.Male.chrY.RData ]; then
    rm ${untransformed_methylation_adjusted}.Male.chrY.[0-9]*.RData
fi

echo "Successfully aggregating methylation_adjustment2 chunks"

