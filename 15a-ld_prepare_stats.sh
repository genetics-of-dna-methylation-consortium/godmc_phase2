#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_15a_logfile})
print_version

mkdir -p ${ld_dir}
mkdir -p ${ld_blocks_dir}
mkdir -p ${ld_prepare_dir}
mkdir -p ${ld_hail_tmp_dir}

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

hail_runtime_args=""
if [ -n "${ld_hail_local_cores}" ]; then
	hail_runtime_args="${hail_runtime_args} --hail-local-cores ${ld_hail_local_cores}"
fi
if [ -n "${ld_hail_driver_memory_gb}" ]; then
	hail_runtime_args="${hail_runtime_args} --hail-driver-memory-gb ${ld_hail_driver_memory_gb}"
fi

python resources/genetics/ld_prepare_stats.py \
	--study-name ${study_name} \
	--bfile ${bfile} \
	--covariates ${covariates_intersect} \
	--output-dir ${ld_prepare_dir} \
	--log-file ${section_15a_logfile} \
	--chromosome ${ld_chromosome} \
	--hail-partitions ${ld_hail_partitions} \
	--hail-tmp-dir ${ld_hail_tmp_dir} \
	--a-block-size ${ld_a_block_size} \
	--a-chunk-rows ${ld_a_chunk_rows} \
	--a-max-dense-gb ${ld_a_max_dense_gb} \
	${hail_runtime_args}

echo "Successfully prepared LD cohort scaffold outputs"
