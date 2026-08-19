#!/usr/bin/env bash
# ld_pack.sh — the single source of truth for section-15 (LD) archive packing.
#
# Source this file and call ld_pack_archive to turn one or more files/dirs into
# the tar+md5+GPG artifact triple that the section-15 upload path ships and that
# ld_reassemble_cohort.py consumes. Both 15c-ld_run_upload.sh (streaming,
# per-chunk) and ld_encrypt_cohort.sh (whole per-chromosome dir) drive it, so the
# tar/md5/encrypt convention lives in exactly one place.
#
# ld_pack_archive <out_dir> <base> <src_parent> <member>...
#   Produces, under <out_dir>:
#     <base>.tgz.aes   AES256 symmetric GPG encryption of the gzipped tar
#     <base>.md5sum    md5 of the *plaintext* tar (verified on the decrypt side)
#   <member>... are paths relative to <src_parent> (one for a chunk, several for
#   a scaffold). The intermediate <base>.tgz is removed once encrypted.
#
#   Resume guard: if both <base>.tgz.aes and <base>.md5sum already exist the call
#   is a no-op (prints "[ld_pack] skip <base>"), so a re-run after a partial
#   failure does not re-encrypt already-staged artifacts.
#
#   Set LD_GPG_PASSPHRASE_FILE or ld_gpg_passphrase_file to a readable file
#   containing the shared encryption passphrase. This is required so section-15
#   packaging can run safely under batch schedulers without pinentry.
#
#   Override the gpg binary/wrapper for testing via the GPG env var.

ld_gpg_passphrase_file_path () {
	if [ -n "${LD_GPG_PASSPHRASE_FILE:-}" ]; then
		echo "${LD_GPG_PASSPHRASE_FILE}"
	else
		echo "${ld_gpg_passphrase_file:-}"
	fi
}

ld_require_gpg_passphrase_file () {
	local passphrase_file="$1"
	if [ -z "${passphrase_file}" ]; then
		echo "Problem: section-15 GPG packaging requires a passphrase file." >&2
		echo "Set ld_gpg_passphrase_file in config, or export LD_GPG_PASSPHRASE_FILE." >&2
		echo "Example: ld_gpg_passphrase_file=\"/secure/path/godmc15.pass\"" >&2
		return 1
	fi
	if [ ! -r "${passphrase_file}" ]; then
		echo "Problem: LD GPG passphrase file is not readable: ${passphrase_file}" >&2
		return 1
	fi
	if [ ! -s "${passphrase_file}" ]; then
		echo "Problem: LD GPG passphrase file is empty: ${passphrase_file}" >&2
		return 1
	fi
}

ld_verify_archive () {
	local out_dir="$1" base="$2"
	local gpg="${GPG:-gpg}" passphrase_file expected actual
	passphrase_file="$(ld_gpg_passphrase_file_path)"
	ld_require_gpg_passphrase_file "${passphrase_file}" || return 1
	if [ ! -f "${out_dir}/${base}.tgz.aes" ] || [ ! -f "${out_dir}/${base}.md5sum" ]; then
		echo "Problem: incomplete encrypted archive pair for ${base}" >&2
		return 1
	fi
	expected="$(cut -d ' ' -f 1 "${out_dir}/${base}.md5sum")"
	if ! actual="$(
		set -o pipefail
		"${gpg}" --batch --yes --pinentry-mode loopback \
			--passphrase-file "${passphrase_file}" \
			--decrypt "${out_dir}/${base}.tgz.aes" 2>/dev/null | md5sum | cut -d ' ' -f 1
	)"; then
		echo "Problem: GPG decryption failed while verifying ${out_dir}/${base}.tgz.aes" >&2
		return 1
	fi
	if [ -z "${expected}" ] || [ "${actual}" != "${expected}" ]; then
		echo "Problem: encrypted archive verification failed for ${out_dir}/${base}.tgz.aes" >&2
		return 1
	fi
}

ld_pack_archive () {
	local out_dir="$1" base="$2" src_parent="$3"
	shift 3
	local gpg="${GPG:-gpg}"
	local passphrase_file member source_available=1

	mkdir -p "${out_dir}"
	if [ -f "${out_dir}/${base}.tgz.aes" ] && [ -f "${out_dir}/${base}.md5sum" ]; then
		if ld_verify_archive "${out_dir}" "${base}"; then
			echo "[ld_pack] skip ${base} (already packed and verified)"
			return 0
		fi
		for member in "$@"; do
			if [ ! -e "${src_parent}/${member}" ]; then
				source_available=0
			fi
		done
		if [ "${source_available}" -ne 1 ]; then
			echo "Problem: existing encrypted archive is invalid and source data is unavailable for ${base}" >&2
			return 1
		fi
		echo "[ld_pack] replacing invalid existing archive ${base}" >&2
		rm -f "${out_dir}/${base}.tgz.aes" "${out_dir}/${base}.md5sum"
	fi
	passphrase_file="$(ld_gpg_passphrase_file_path)"
	ld_require_gpg_passphrase_file "${passphrase_file}" || return 1

	rm -f "${out_dir}/.uploaded_${base}.tgz.aes" \
		"${out_dir}/.uploaded_${base}.md5sum"
	if ! tar czf "${out_dir}/${base}.tgz" -C "${src_parent}" "$@"; then
		rm -f "${out_dir}/${base}.tgz"
		echo "Problem: failed to create archive ${base}.tgz" >&2
		return 1
	fi
	if ! ( cd "${out_dir}" && md5sum "${base}.tgz" > "${base}.md5sum" ); then
		rm -f "${out_dir}/${base}.tgz" "${out_dir}/${base}.md5sum"
		echo "Problem: failed to checksum archive ${base}.tgz" >&2
		return 1
	fi
	if ! "${gpg}" --batch --yes --pinentry-mode loopback \
		--passphrase-file "${passphrase_file}" \
		--output "${out_dir}/${base}.tgz.aes" \
		--symmetric --cipher-algo AES256 "${out_dir}/${base}.tgz"; then
		rm -f "${out_dir}/${base}.tgz.aes" \
			"${out_dir}/${base}.tgz" \
			"${out_dir}/${base}.md5sum"
		echo "Problem: GPG encryption failed for ${out_dir}/${base}.tgz" >&2
		echo "If GPG reports loopback pinentry is not allowed, add 'allow-loopback-pinentry' to ~/.gnupg/gpg-agent.conf and restart gpg-agent." >&2
		return 1
	fi
	if ! ld_verify_archive "${out_dir}" "${base}"; then
		rm -f "${out_dir}/${base}.tgz.aes" \
			"${out_dir}/${base}.md5sum" \
			"${out_dir}/${base}.tgz"
		echo "Problem: GPG archive round-trip verification failed for ${base}" >&2
		return 1
	fi
	rm -f "${out_dir}/${base}.tgz"
}
