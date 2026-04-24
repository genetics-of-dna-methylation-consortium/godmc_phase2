#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_15b_logfile})
print_version

mkdir -p ${ld_dir}
mkdir -p ${ld_blocks_dir}
mkdir -p ${ld_aggregate_dir}

if [ ! -f "${ld_prepare_dir}/manifest.json" ]
then
	echo "Problem: cohort LD scaffold manifest is required at ${ld_prepare_dir}/manifest.json"
	exit 1
fi

python resources/genetics/ld_aggregate_stats.py \
	--study-name ${study_name} \
	--cohort-dir ${ld_prepare_dir} \
	--output-dir ${ld_aggregate_dir} \
	--log-file ${section_15b_logfile}

echo "Successfully prepared LD aggregation scaffold outputs"
