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
	mkdir -p "${out_dir}"
	tar czf "${out_dir}/${base}.tgz" -C "${src_parent}" "${member}"
	( cd "${out_dir}" && md5sum "${base}.tgz" > "${base}.md5sum" )
	"${GPG}" --output "${out_dir}/${base}.tgz.aes" \
		--symmetric --cipher-algo AES256 "${out_dir}/${base}.tgz"
	rm -f "${out_dir}/${base}.tgz"
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

main () {
	source resources/setup.sh "$@"
	set -- $concatenated
	echo "[15c] main is implemented in a later task"
}

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
	main "$@"
fi
