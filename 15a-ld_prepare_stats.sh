#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_15a_logfile})
print_version

mkdir -p ${ld_dir}
mkdir -p ${ld_blocks_dir}
mkdir -p ${ld_prepare_dir}

if [ ! -f "${bfile}.bed" ] || [ ! -f "${bfile}.bim" ] || [ ! -f "${bfile}.fam" ]
then
	echo "Problem: cleaned section-02 genotype files are required at ${bfile}"
	exit 1
fi

if [ ! -f "${covariates_intersect}" ]
then
	echo "Problem: intersected covariates file is required at ${covariates_intersect}"
	exit 1
fi

python resources/genetics/ld_prepare_stats.py \
	--study-name ${study_name} \
	--bfile ${bfile} \
	--covariates ${covariates_intersect} \
	--output-dir ${ld_prepare_dir} \
	--log-file ${section_15a_logfile}

echo "Successfully prepared LD cohort scaffold outputs"
