#!/usr/bin/env bash

ld_validate_chromosome_15 () {
	local chr="$1"
	if ! [[ "${chr}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
		echo "Problem: invalid section-15 chromosome '${chr}'; expected 1-22" >&2
		return 1
	fi
}

ld_target_chromosomes_15 () {
	ld_validate_chromosome_15 "${ld_chromosome}" || return 1
	echo "${ld_chromosome}"
}

ld_resolve_chromosome_dir_15 () {
	local base_dir="$1" chr="$2"
	if [ "$(basename "${base_dir}")" = "chr${chr}" ]; then
		echo "${base_dir}"
	else
		echo "${base_dir}/chr${chr}"
	fi
}

ld_upload_receipt_present_15 () {
	local out_dir="$1" base="$2"
	[ -f "${out_dir}/.uploaded_${base}.tgz.aes" ] && \
		[ -f "${out_dir}/.uploaded_${base}.md5sum" ]
}

ld_upload_artefact_ready_15 () {
	local out_dir="$1" base="$2"
	if ld_upload_receipt_present_15 "${out_dir}" "${base}"; then
		return 0
	fi
	[ -f "${out_dir}/${base}.tgz.aes" ] && [ -f "${out_dir}/${base}.md5sum" ]
}
