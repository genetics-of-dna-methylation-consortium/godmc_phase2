#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated
set -o pipefail
source "${scripts_directory}/resources/genetics/ld_workflow.sh"

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
		${hail_runtime_args} || return 1

	check_prepared_outputs "${outdir}" "${chr}" || return 1
	touch "${outdir}/.prepared" || return 1

	# Hail materialises the chromosome-wide genotype BlockMatrix under the temp
	# directory; it is not reused across chromosomes, so reclaim it now to keep
	# peak local disk to one chromosome.
	if [ -n "${ld_hail_tmp_dir}" ] && [ -d "${ld_hail_tmp_dir}" ]; then
		rm -rf "${ld_hail_tmp_dir}" || return 1
	fi

	echo "Successfully prepared LD cohort scaffold outputs for chr${chr}"
}

run_chromosome () {
	local chr="$1" outdir="$2" logfile="$3"
	if [ -f "${outdir}/.prepared" ] && check_prepared_outputs "${outdir}" "${chr}"; then
		echo "[15a] chr${chr} already prepared; skip"
		echo "Successfully prepared LD cohort scaffold outputs for chr${chr}"
		return 0
	fi
	run_prepare_stats "${chr}" "${outdir}" "${logfile}" || return 1
}

check_inputs

mapfile -t target_chromosomes < <(ld_target_chromosomes_15)
if [ "${#target_chromosomes[@]}" -eq 0 ]; then
	echo "Problem: no section-15 chromosomes selected" >&2
	exit 1
fi
for chr in "${target_chromosomes[@]}"; do
	chromosome_dir="$(ld_resolve_chromosome_dir_15 "${ld_prepare_dir}" "${chr}")"
	chromosome_log="${section_15_dir}/logs_a/chr${chr}.log"
	mkdir -p "$(dirname "${chromosome_log}")"
	if ! run_chromosome "${chr}" "${chromosome_dir}" "${chromosome_log}" 2>&1 | tee "${chromosome_log}"; then
		exit 1
	fi
done
