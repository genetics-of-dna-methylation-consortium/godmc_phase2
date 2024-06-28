#!/bin/bash

set -e

# Initialize variables
config_file="./config"

# Parse options using getopts
while getopts "c:" opt; do
    case $opt in
        c) config_file=$OPTARG ;;
        *) echo "Usage: $0 -c <config_file>"
           exit 1 ;;
    esac
done

# Shift option arguments, so $1 becomes the first positional argument
shift $((OPTIND - 1))

set -e
echo "-----------------------------------------------"
echo ""
echo "Using config located at:" ${config_file}
echo ""
echo "-----------------------------------------------"
	
source ${config_file}
batch_number=${1}
exec &> >(tee ${section_03e_logfile}${batch_number})
print_version

echo "Adjusting methylation for meth PCs"

if [ -n "${1}" ]
then
	re='^[0-9]+$'
	if ! [[ $batch_number =~ $re ]] ; then
		echo "error: Batch variable is not a number"
		exit 1
	fi
	i=${1}
	echo "Running batch ${i} of ${meth_chunks}"
else
	i="NA"
	echo "Running entire set on a single node using ${nthreads} threads."
fi

${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${transformed_methylation_adjusted}.RData \
	${nongenetic_meth_pcs_transformed}.txt \
	${transformed_methylation_adjusted_pcs} \
	transformed \
	${nthreads} \
	${meth_chunks} \
	${i}

if [ ${batch_number} == 1 ] && [ -f ${transformed_methylation_adjusted}.Female.chrX.RData ]; then
${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${transformed_methylation_adjusted}.Female.chrX.RData \
	${nongenetic_meth_pcs_transformed}.txt \
	${transformed_methylation_adjusted_pcs}.Female.chrX \
	transformed \
	${nthreads} \
	${meth_chunks}
fi

if [ ${batch_number} == 1 ] && [ -f ${transformed_methylation_adjusted}.Male.chrX.RData ]; then
${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${transformed_methylation_adjusted}.Male.chrX.RData \
	${nongenetic_meth_pcs_transformed}.txt \
	${transformed_methylation_adjusted_pcs}.Male.chrX \
	transformed \
	${nthreads} \
	${meth_chunks}
fi

if [ ${batch_number} == 1 ] && [ -f ${transformed_methylation_adjusted}.Male.chrY.RData ]; then
${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${transformed_methylation_adjusted}.Male.chrY.RData \
	${nongenetic_meth_pcs_transformed}.txt \
	${transformed_methylation_adjusted_pcs}.Male.chrY \
	transformed \
	${nthreads} \
	${meth_chunks}
fi

${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${untransformed_methylation_adjusted}.RData \
	${nongenetic_meth_pcs_untransformed}.txt \
	${untransformed_methylation_adjusted_pcs} \
	untransformed \
	${nthreads} \
	${meth_chunks} \
	${i}

if [ ${batch_number} == 1 ] && [ -f ${untransformed_methylation_adjusted}.Female.chrX.RData ]; then
${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${untransformed_methylation_adjusted}.Female.chrX.RData \
	${nongenetic_meth_pcs_untransformed}.txt \
	${untransformed_methylation_adjusted_pcs}.Female.chrX \
	untransformed \
	${nthreads} \
	${meth_chunks}
fi

if [ ${batch_number} == 1 ] && [ -f ${untransformed_methylation_adjusted}.Male.chrX.RData ]; then
${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${untransformed_methylation_adjusted}.Male.chrX.RData \
	${nongenetic_meth_pcs_untransformed}.txt \
	${untransformed_methylation_adjusted_pcs}.Male.chrX \
	untransformed \
	${nthreads} \
	${meth_chunks}
fi

if [ ${batch_number} == 1 ] && [ -f ${untransformed_methylation_adjusted}.Male.chrY.RData ]; then
${R_directory}Rscript resources/methylation/adjust_nongeneticpcs.R \
	${untransformed_methylation_adjusted}.Male.chrY.RData \
	${nongenetic_meth_pcs_untransformed}.txt \
	${untransformed_methylation_adjusted_pcs}.Male.chrY \
	untransformed \
	${nthreads} \
	${meth_chunks}
fi

echo "Successfully adjusting non-genetic PCs"
