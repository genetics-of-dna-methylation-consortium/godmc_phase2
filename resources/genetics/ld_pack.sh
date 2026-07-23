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
#   Override the gpg binary/wrapper for testing via the GPG env var.

ld_pack_archive () {
	local out_dir="$1" base="$2" src_parent="$3"
	shift 3
	local gpg="${GPG:-gpg}"

	mkdir -p "${out_dir}"
	if [ -f "${out_dir}/${base}.tgz.aes" ] && [ -f "${out_dir}/${base}.md5sum" ]; then
		echo "[ld_pack] skip ${base} (already packed)"
		return 0
	fi

	tar czf "${out_dir}/${base}.tgz" -C "${src_parent}" "$@"
	( cd "${out_dir}" && md5sum "${base}.tgz" > "${base}.md5sum" )
	"${gpg}" --output "${out_dir}/${base}.tgz.aes" \
		--symmetric --cipher-algo AES256 "${out_dir}/${base}.tgz"
	rm -f "${out_dir}/${base}.tgz"
}
