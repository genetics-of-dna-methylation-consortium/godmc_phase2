#!/bin/bash

set -e
source config
exec &> >(tee ${section_03b_logfile}_${1}_aggregation)

subset=$1

echo "Aggregating methylation_adjustment1 chunks"
${R_directory}Rscript ${scripts_directory}/resources/methylation/aggregate_chunks.R \
	${transformed_methylation_adjusted} \
	${meth_chunks} \
	${methylation_array} \
	${methylation_no_outliers} \
	${subset}

${R_directory}Rscript ${scripts_directory}/resources/methylation/aggregate_chunks.R \
	${untransformed_methylation_adjusted} \
	${meth_chunks} \
        ${methylation_array} \
        ${methylation_no_outliers} \
	${subset}

if [ "${related}" = "yes" ]
then
    echo "You have specified that the data is family data. Aggregating the results from polygenic effects"
    ${R_directory}Rscript ${scripts_directory}/resources/methylation/aggregate_chunk_classes.R \
   	${transformed_methylation_adjusted} \
	${meth_chunks} \
	${section_03_dir}/${subset}_transformed_classes.RData \
	${subset}

    ${R_directory}Rscript ${scripts_directory}/resources/methylation/aggregate_chunk_classes.R \
        ${untransformed_methylation_adjusted} \
        ${meth_chunks} \
        ${section_03_dir}/${subset}_untransformed_classes.RData \
	${subset}
fi

echo "Succesfully aggregating methylation_adjustment1 chunks"
