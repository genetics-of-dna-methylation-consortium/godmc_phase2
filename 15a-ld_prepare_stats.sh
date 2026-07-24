#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated
set -o pipefail

mkdir -p "${ld_dir}" "${ld_blocks_dir}" "${ld_prepare_dir}" "${ld_hail_tmp_dir}"

check_inputs () {
	if [ ! -f "${bfile}.bed" ] || [ ! -f "${bfile}.bim" ] || [ ! -f "${bfile}.fam" ]
	then
		echo "Problem: cleaned section-02 genotype files are required at ${bfile}"
		exit 1
	fi

	if [ ! -f "${covariates_intersect}" ]
	then
		echo "Problem: section-03a mQTL-aligned covariates file is required at ${covariates_intersect}"
		echo "Please run 03a-methylation_variables.sh before 15a-ld_prepare_stats.sh"
		exit 1
	fi
}

check_prepared_outputs () {
	local outdir="$1" chr="$2" f
	for f in manifest.json variants.tsv.gz B.npy D.npy checksums.json qc_report.txt; do
		if [ ! -f "${outdir}/${f}" ]; then
			return 1
		fi
	done
	[ -d "${outdir}/A_blocks/chr${chr}" ]
}

resolve_chromosome_dir () {
	local chr="$1"
	if [ "$(basename "${ld_prepare_dir}")" = "chr${chr}" ]; then
		echo "${ld_prepare_dir}"
	else
		echo "${ld_prepare_dir}/chr${chr}"
	fi
}

run_prepare_stats () {
	local chr="$1" outdir="$2" logfile="$3"
	local hail_runtime_args=""
	mkdir -p "${outdir}" "$(dirname "${logfile}")"

	if [ -n "${ld_hail_local_cores}" ]; then
		hail_runtime_args="${hail_runtime_args} --hail-local-cores ${ld_hail_local_cores}"
	fi
	if [ -n "${ld_hail_driver_memory_gb}" ]; then
		hail_runtime_args="${hail_runtime_args} --hail-driver-memory-gb ${ld_hail_driver_memory_gb}"
	fi

	print_version
	python resources/genetics/ld_prepare_stats.py \
		--study-name "${study_name}" \
		--bfile "${bfile}" \
		--covariates "${covariates_intersect}" \
		--output-dir "${outdir}" \
		--log-file "${logfile}" \
		--chromosome "${chr}" \
		--hail-partitions "${ld_hail_partitions}" \
		--hail-tmp-dir "${ld_hail_tmp_dir}" \
		--a-block-size "${ld_a_block_size}" \
		--a-chunk-rows "${ld_a_chunk_rows}" \
		--a-max-dense-gb "${ld_a_max_dense_gb}" \
		${hail_runtime_args}

	check_prepared_outputs "${outdir}" "${chr}"
	touch "${outdir}/.prepared"
	echo "Successfully prepared LD cohort scaffold outputs for chr${chr}"
}

run_chromosome () {
	local chr="$1" outdir="$2" logfile="$3"
	if [ -f "${outdir}/.prepared" ] && check_prepared_outputs "${outdir}" "${chr}"; then
		echo "[15a] chr${chr} already prepared; skip"
		echo "Successfully prepared LD cohort scaffold outputs for chr${chr}"
		return 0
	fi
	run_prepare_stats "${chr}" "${outdir}" "${logfile}"
}

if [ "${ld_chromosome}" = "all" ]; then
	echo "Problem: section 15 must be prepared one chromosome at a time to keep disk use bounded."
	echo "Run, for example: ld_chromosome=22 bash 15a-ld_prepare_stats.sh -c config"
	exit 1
fi

check_inputs

chromosome_dir="$(resolve_chromosome_dir "${ld_chromosome}")"
chromosome_log="${section_15_dir}/logs_a/chr${ld_chromosome}.log"
mkdir -p "$(dirname "${chromosome_log}")"
exec &> >(tee "${chromosome_log}")

run_chromosome "${ld_chromosome}" "${chromosome_dir}" "${chromosome_log}"
