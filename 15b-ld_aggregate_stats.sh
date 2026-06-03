#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_15b_logfile})
print_version

mkdir -p ${ld_dir}
mkdir -p ${ld_precursor_dir}
mkdir -p ${ld_panel_dir}

if [ "${ld_15b_mode}" = "accumulate" ]; then
	if [ ! -f "${ld_prepare_dir}/manifest.json" ]; then
		echo "Problem: cohort 15a manifest is required at ${ld_prepare_dir}/manifest.json"
		exit 1
	fi
	python resources/genetics/ld_aggregate_stats.py \
		--mode accumulate \
		--cohort-dir ${ld_prepare_dir} \
		--precursor-dir ${ld_precursor_dir} \
		--log-file ${section_15b_logfile}
	echo "Successfully accumulated cohort into the LD precursor"
elif [ "${ld_15b_mode}" = "finalise" ]; then
	min_cohorts_arg=""
	if [ -n "${ld_min_cohorts}" ]; then
		min_cohorts_arg="--min-cohorts ${ld_min_cohorts}"
	fi
	python resources/genetics/ld_aggregate_stats.py \
		--mode finalise \
		--precursor-dir ${ld_precursor_dir} \
		--panel-dir ${ld_panel_dir} \
		--maf-threshold ${ld_maf_threshold} \
		--min-adj-diag ${ld_min_adj_diag} \
		${min_cohorts_arg} \
		--log-file ${section_15b_logfile}
	echo "Successfully finalised the pooled LD panel"
else
	echo "Problem: ld_15b_mode must be 'accumulate' or 'finalise'; got '${ld_15b_mode}'"
	exit 1
fi
