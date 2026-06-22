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

main () {
	source resources/setup.sh "$@"
	set -- $concatenated
	echo "[15c] main is implemented in a later task"
}

if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
	main "$@"
fi
