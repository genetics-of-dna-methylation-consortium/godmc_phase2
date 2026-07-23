#!/usr/bin/env bash
# 15c-ld_run_upload.sh — cohort-side section-15 (LD) run -> encrypt -> upload
# (Imperial) -> verify -> delete, one chromosome at a time, so peak disk stays
# at roughly one chromosome's A_blocks (~100 GB) instead of the ~1.3 TB a
# whole-genome run would require.
#
# Usage: ./15c-ld_run_upload.sh -c <config>
#
# Resumable: each chromosome writes results/15/cohort_stats/chr<C>/.uploaded
# only after every artifact has shipped; re-running skips completed chromosomes.
#
# Test/override hooks (no effect in normal operation):
#   GPG            gpg binary/wrapper (default: gpg)
#   LD_SHIP_CMD    command taking one file path; replaces the rsync transport
#   LD_PREPARE_CMD command "<chr> <outdir>"; replaces the 15a compute step
set -euo pipefail

GPG="${GPG:-gpg}"

# Shared tar+md5+GPG primitive (single source of truth for section-15 packing).
source "$(dirname "${BASH_SOURCE[0]}")/resources/genetics/ld_pack.sh"

# chr_done <outdir> — true if the chromosome's upload sentinel exists.
chr_done () {
	[ -f "$1/.uploaded" ]
}

# ship_file <path> — transfer one file to Imperial (or via the test override).
ship_file () {
	local path="$1"
	if [ -n "${LD_SHIP_CMD:-}" ]; then
		"${LD_SHIP_CMD}" "${path}"
		return $?
	fi
	rsync --partial --append --checksum \
		-e "ssh -i ${imperial_key}" \
		"${path}" "${imperial_user}@${imperial_host}:${imperial_path}/"
}

# encrypt_archive <src_parent> <member> <out_dir> <base>
# tar+md5(plaintext)+gpg one path; removes the intermediate .tgz.
encrypt_archive () {
	local src_parent="$1" member="$2" out_dir="$3" base="$4"
	ld_pack_archive "${out_dir}" "${base}" "${src_parent}" "${member}"
}

# process_chunk <chunk_dir> <ablocks_dir> <out_dir> <study_tag>
# encrypt -> ship .aes + .md5sum -> on success delete source + staged artifacts;
# on ship failure leave everything in place and return 1.
process_chunk () {
	local chunk_dir="$1" ablocks_dir="$2" out_dir="$3" study_tag="$4"
	local chr_name chunk_name base
	chr_name="$(basename "$(dirname "${chunk_dir}")")"
	chunk_name="$(basename "${chunk_dir}")"
	base="${study_tag}_15_${chr_name}_${chunk_name}"
	encrypt_archive "${ablocks_dir}" "${chr_name}/${chunk_name}" "${out_dir}" "${base}"
	if ship_file "${out_dir}/${base}.tgz.aes" && ship_file "${out_dir}/${base}.md5sum"; then
		rm -rf "${chunk_dir}" "${out_dir}/${base}.tgz.aes" "${out_dir}/${base}.md5sum"
		echo "[15c] shipped ${base}"
		return 0
	fi
	echo "[15c] ship FAILED for ${base}; leaving source in place" >&2
	return 1
}

# prepare_chromosome <C> <outdir> — run 15a compute for chromosome C into outdir.
prepare_chromosome () {
	local chr="$1" outdir="$2"
	mkdir -p "${outdir}"
	if [ -n "${LD_PREPARE_CMD:-}" ]; then
		"${LD_PREPARE_CMD}" "${chr}" "${outdir}"
		return $?
	fi
	local hail_runtime_args=""
	if [ -n "${ld_hail_local_cores:-}" ]; then
		hail_runtime_args="${hail_runtime_args} --hail-local-cores ${ld_hail_local_cores}"
	fi
	if [ -n "${ld_hail_driver_memory_gb:-}" ]; then
		hail_runtime_args="${hail_runtime_args} --hail-driver-memory-gb ${ld_hail_driver_memory_gb}"
	fi
	python "${scripts_directory}/resources/genetics/ld_prepare_stats.py" \
		--study-name "${study_name}" \
		--bfile "${bfile}" \
		--covariates "${covariates_intersect}" \
		--output-dir "${outdir}" \
		--log-file "${section_15a_logfile}" \
		--chromosome "${chr}" \
		--hail-partitions "${ld_hail_partitions}" \
		--hail-tmp-dir "${ld_hail_tmp_dir}" \
		--a-block-size "${ld_a_block_size}" \
		--a-chunk-rows "${ld_a_chunk_rows}" \
		--a-max-dense-gb "${ld_a_max_dense_gb}" \
		${hail_runtime_args}
}

# check_chromosome <outdir> — required scaffold + A_blocks present.
check_chromosome () {
	local outdir="$1" f
	for f in manifest.json variants.tsv.gz D.npy B.npy; do
		if [ ! -f "${outdir}/${f}" ]; then
			echo "[15c] check failed: missing ${f} in ${outdir}" >&2
			return 1
		fi
	done
	if [ ! -d "${outdir}/A_blocks" ]; then
		echo "[15c] check failed: missing A_blocks in ${outdir}" >&2
		return 1
	fi
}

# ship_scaffold <outdir> <out_dir> <study_tag>
ship_scaffold () {
	local outdir="$1" out_dir="$2" study_tag="$3"
	local base="${study_tag}_15_scaffold"
	ld_pack_archive "${out_dir}" "${base}" "${outdir}" \
		manifest.json variants.tsv.gz D.npy B.npy checksums.json qc_report.txt
	if ship_file "${out_dir}/${base}.tgz.aes" && ship_file "${out_dir}/${base}.md5sum"; then
		rm -f "${out_dir}/${base}.tgz.aes" "${out_dir}/${base}.md5sum"
		echo "[15c] shipped ${base}"
		return 0
	fi
	echo "[15c] ship FAILED for ${base}; leaving source in place" >&2
	return 1
}

# process_chromosome <C> <study> <cohort_root> <out_dir>
process_chromosome () {
	local chr="$1" study="$2" cohort_root="$3" out_dir="$4"
	local outdir="${cohort_root}/chr${chr}"
	if chr_done "${outdir}"; then
		echo "[15c] chr${chr} already uploaded; skip"
		return 0
	fi
	prepare_chromosome "${chr}" "${outdir}"
	check_chromosome "${outdir}" || return 1
	local ablocks="${outdir}/A_blocks" tag="${study}_chr${chr}"
	local chunk_dir
	shopt -s nullglob
	for chunk_dir in "${ablocks}"/chr*/chunk_*; do
		[ -d "${chunk_dir}" ] || continue
		if ! process_chunk "${chunk_dir}" "${ablocks}" "${out_dir}" "${tag}"; then
			shopt -u nullglob
			return 1
		fi
	done
	shopt -u nullglob
	ship_scaffold "${outdir}" "${out_dir}" "${tag}" || return 1
	touch "${outdir}/.uploaded"
	rm -rf "${ablocks}"
	echo "[15c] chr${chr} complete"
}

main () {
	source resources/setup.sh "$@"
	set -- $concatenated

	exec &> >(tee "${section_15a_logfile}")
	print_version

	local cohort_root="${ld_prepare_dir}"
	local upload_dir="${section_15_dir}/upload"
	mkdir -p "${cohort_root}" "${upload_dir}" "${ld_hail_tmp_dir}"

	if [ ! -f "${bfile}.bed" ]; then
		echo "Problem: cleaned section-02 genotype files are required at ${bfile}"
		exit 1
	fi
	if [ ! -f "${covariates_intersect}" ]; then
		echo "Problem: section-03a mQTL-aligned covariates are required at ${covariates_intersect}"
		exit 1
	fi

	local failed=0 chr
	for chr in ${ld_chromosomes}; do
		echo "[15c] ===== chromosome ${chr} ====="
		if ! process_chromosome "${chr}" "${study_name}" "${cohort_root}" "${upload_dir}"; then
			echo "[15c] chromosome ${chr} did NOT complete; re-run to resume" >&2
			failed=1
		fi
	done

	if [ "${failed}" -ne 0 ]; then
		echo "[15c] one or more chromosomes failed; re-run ./15c-ld_run_upload.sh to resume"
		exit 1
	fi
	echo "Successfully ran and uploaded all section-15 LD cohort chromosomes"
}

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
	main "$@"
fi
