#!/usr/bin/env bash

source resources/setup.sh "$@"
set -- $concatenated
set -o pipefail

source "${scripts_directory}/resources/genetics/ld_pack.sh"
source "${scripts_directory}/resources/genetics/ld_workflow.sh"

section_15_upload_dir="${section_15_dir}/upload"

mkdir -p "${ld_prepare_dir}" "${section_15_upload_dir}" "${section_15_dir}/logs_b"

required_files="manifest.json variants.tsv.gz B.npy D.npy checksums.json qc_report.txt"

check_required_files () {
	local outdir="$1" chr="$2" f
	for f in ${required_files}; do
		if [ ! -f "${outdir}/${f}" ]; then
			echo "Problem: missing ${f} in ${outdir}" >&2
			return 1
		fi
	done
	if [ ! -d "${outdir}/A_blocks/chr${chr}" ] && [ ! -f "${outdir}/.packaged" ]; then
		echo "Problem: missing ${outdir}/A_blocks/chr${chr}" >&2
		return 1
	fi
}

manifest_chunks () {
	local outdir="$1" chr="$2"
	PYTHONPATH="${scripts_directory}/resources/genetics:${PYTHONPATH:-}" python - "${outdir}" "${chr}" <<'PY'
import json
import sys
from pathlib import Path

outdir = Path(sys.argv[1])
chrom = sys.argv[2]
manifest = json.loads((outdir / "manifest.json").read_text(encoding="utf-8"))
actual_filter = manifest.get("variant_index", {}).get("chromosome_filter")
if str(actual_filter) != str(chrom):
    raise SystemExit(
        f"manifest chromosome_filter is {actual_filter!r}, expected {chrom!r}"
    )
chrom_meta = manifest.get("A_blocks", {}).get("chromosomes", {}).get(str(chrom))
if not chrom_meta:
    raise SystemExit(f"manifest lacks A_blocks metadata for chr{chrom}")
chunks = [chunk.get("name") for chunk in chrom_meta.get("chunks", [])]
if any(not name for name in chunks):
    raise SystemExit(f"manifest has an unnamed A-block chunk for chr{chrom}")
if int(chrom_meta.get("n_chunks", -1)) != len(chunks):
    raise SystemExit("manifest n_chunks does not match chunk list")
if not chunks:
    raise SystemExit(f"manifest lists no chunks for chr{chrom}")
for name in chunks:
    print(name)
PY
}

verify_checksum_members () {
	local outdir="$1"
	shift
	PYTHONPATH="${scripts_directory}/resources/genetics:${PYTHONPATH:-}" python - "${outdir}" "$@" <<'PY'
import sys
from pathlib import Path
import ld_checksums

outdir = Path(sys.argv[1])
prefixes = sys.argv[2:]
document = ld_checksums.read_cohort_checksums(outdir)
files = document.get("files") or {}
selected = {
    rel: digest
    for rel, digest in files.items()
    if any(rel == prefix or rel.startswith(prefix.rstrip("/") + "/") for prefix in prefixes)
}
if not selected:
    raise SystemExit(f"checksums.json contains no entries for {', '.join(prefixes)}")
for rel, expected in sorted(selected.items()):
    path = outdir / rel
    if not path.is_file():
        raise SystemExit(f"checksummed artefact is missing: {rel}")
    actual = ld_checksums.hash_file(path)
    if actual != expected:
        raise SystemExit(f"checksum mismatch for {rel}: expected {expected}, got {actual}")
PY
}

artifact_pair_exists () {
	local base="$1"
	[ -f "${section_15_upload_dir}/${base}.tgz.aes" ] && [ -f "${section_15_upload_dir}/${base}.md5sum" ]
}

artifact_pair_valid () {
	local base="$1"
	artifact_pair_exists "${base}" && ld_verify_archive "${section_15_upload_dir}" "${base}"
}

all_artifacts_exist () {
	local chr="$1" outdir="$2" chunk base chunks
	base="${study_name}_chr${chr}_15_scaffold"
	artifact_pair_exists "${base}" || return 1
	chunks="$(manifest_chunks "${outdir}" "${chr}")" || return 1
	[ -n "${chunks}" ] || return 1
	while IFS= read -r chunk; do
		base="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
		artifact_pair_exists "${base}" || return 1
	done <<< "${chunks}"
}

all_artifacts_valid () {
	local chr="$1" outdir="$2" chunk base chunks
	base="${study_name}_chr${chr}_15_scaffold"
	artifact_pair_valid "${base}" || return 1
	chunks="$(manifest_chunks "${outdir}" "${chr}")" || return 1
	[ -n "${chunks}" ] || return 1
	while IFS= read -r chunk; do
		base="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
		artifact_pair_valid "${base}" || return 1
	done <<< "${chunks}"
}

pack_scaffold () {
	local chr="$1" outdir="$2" base
	base="${study_name}_chr${chr}_15_scaffold"
	verify_checksum_members "${outdir}" \
		manifest.json variants.tsv.gz D.npy B.npy || return 1
	rm -f "${section_15_upload_dir}/${base}.tgz.aes" \
		"${section_15_upload_dir}/${base}.md5sum" \
		"${section_15_upload_dir}/.uploaded_${base}.tgz.aes" \
		"${section_15_upload_dir}/.uploaded_${base}.md5sum"
	ld_pack_archive "${section_15_upload_dir}" "${base}" "${outdir}" \
		manifest.json variants.tsv.gz D.npy B.npy checksums.json qc_report.txt || return 1
	artifact_pair_exists "${base}" || return 1
}

pack_chunk () {
	local chr="$1" outdir="$2" chunk="$3" base chunk_dir
	base="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
	chunk_dir="${outdir}/A_blocks/chr${chr}/${chunk}"
	if [ ! -d "${chunk_dir}" ]; then
		if artifact_pair_valid "${base}"; then
			return 0
		fi
		echo "Problem: raw chunk ${chunk_dir} is missing and ${base} has not been packaged correctly" >&2
		return 1
	fi
	verify_checksum_members "${outdir}" "A_blocks/chr${chr}/${chunk}" || return 1
	rm -f "${section_15_upload_dir}/${base}.tgz.aes" \
		"${section_15_upload_dir}/${base}.md5sum" \
		"${section_15_upload_dir}/.uploaded_${base}.tgz.aes" \
		"${section_15_upload_dir}/.uploaded_${base}.md5sum"
	ld_pack_archive "${section_15_upload_dir}" "${base}" "${outdir}/A_blocks" "chr${chr}/${chunk}" || return 1
	artifact_pair_exists "${base}" || return 1
	rm -rf "${chunk_dir}" || return 1
}

check_unexpected_chunks () {
	local chr="$1" outdir="$2" expected_file="$3" chunk_dir chunk
	if [ ! -d "${outdir}/A_blocks/chr${chr}" ]; then
		return 0
	fi
	for chunk_dir in "${outdir}/A_blocks/chr${chr}"/chunk_*; do
		[ -e "${chunk_dir}" ] || continue
		chunk="$(basename "${chunk_dir}")"
		if ! grep -qx "${chunk}" "${expected_file}"; then
			echo "Problem: unexpected chunk directory ${chunk_dir}" >&2
			return 1
		fi
	done
}

process_chromosome () {
	local chr="$1" outdir="$2" tmp_chunks=""
	print_version
	echo "[15b] ===== chromosome ${chr} ====="
	check_required_files "${outdir}" "${chr}" || return 1
	if [ -f "${outdir}/.packaged" ] && [ ! -d "${outdir}/A_blocks" ]; then
		all_artifacts_valid "${chr}" "${outdir}" || return 1
		pack_scaffold "${chr}" "${outdir}" || return 1
		echo "[15b] chr${chr} chunks already packaged and verified; scaffold refreshed"
		echo "Successfully packaged LD cohort chromosome chr${chr}"
		return 0
	fi

	tmp_chunks="$(mktemp)" || return 1
	if ! manifest_chunks "${outdir}" "${chr}" > "${tmp_chunks}"; then
		rm -f "${tmp_chunks}"
		return 1
	fi
	if ! check_unexpected_chunks "${chr}" "${outdir}" "${tmp_chunks}"; then
		rm -f "${tmp_chunks}"
		return 1
	fi
	while IFS= read -r chunk; do
		if ! pack_chunk "${chr}" "${outdir}" "${chunk}"; then
			rm -f "${tmp_chunks}"
			return 1
		fi
	done < "${tmp_chunks}"
	rm -f "${tmp_chunks}" || return 1

	pack_scaffold "${chr}" "${outdir}" || return 1
	all_artifacts_exist "${chr}" "${outdir}" || return 1
	rm -rf "${outdir}/A_blocks" || return 1
	touch "${outdir}/.packaged" || return 1
	echo "Successfully packaged LD cohort chromosome chr${chr}"
}

mapfile -t target_chromosomes < <(ld_target_chromosomes_15)
if [ "${#target_chromosomes[@]}" -eq 0 ]; then
	echo "Problem: no section-15 chromosomes selected" >&2
	exit 1
fi
for chr in "${target_chromosomes[@]}"; do
	chromosome_dir="$(ld_resolve_chromosome_dir_15 "${ld_prepare_dir}" "${chr}")"
	chromosome_log="${section_15_dir}/logs_b/chr${chr}.log"
	mkdir -p "$(dirname "${chromosome_log}")"
	if ! process_chromosome "${chr}" "${chromosome_dir}" 2>&1 | tee "${chromosome_log}"; then
		exit 1
	fi
done
