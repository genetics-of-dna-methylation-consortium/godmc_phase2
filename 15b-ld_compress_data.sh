#!/usr/bin/env bash

source resources/setup.sh "$@"
set -- $concatenated
set -o pipefail

source "${scripts_directory}/resources/genetics/ld_pack.sh"

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

resolve_chromosome_dir () {
	local chr="$1"
	if [ "$(basename "${ld_prepare_dir}")" = "chr${chr}" ]; then
		echo "${ld_prepare_dir}"
	else
		echo "${ld_prepare_dir}/chr${chr}"
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
chunks = [name for name in chunks if name]
if int(chrom_meta.get("n_chunks", -1)) != len(chunks):
    raise SystemExit("manifest n_chunks does not match chunk list")
if not chunks:
    raise SystemExit(f"manifest lists no chunks for chr{chrom}")
for name in chunks:
    print(name)
PY
}

verify_checksums () {
	local outdir="$1"
	PYTHONPATH="${scripts_directory}/resources/genetics:${PYTHONPATH:-}" python - "${outdir}" <<'PY'
import sys
from pathlib import Path
import ld_checksums

ld_checksums.verify_cohort_checksums(Path(sys.argv[1]))
PY
}

artifact_pair_exists () {
	local base="$1"
	[ -f "${section_15_upload_dir}/${base}.tgz.aes" ] && [ -f "${section_15_upload_dir}/${base}.md5sum" ]
}

all_artifacts_exist () {
	local chr="$1" outdir="$2" chunk base
	base="${study_name}_chr${chr}_15_scaffold"
	artifact_pair_exists "${base}" || return 1
	while IFS= read -r chunk; do
		base="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
		artifact_pair_exists "${base}" || return 1
	done < <(manifest_chunks "${outdir}" "${chr}")
}

pack_scaffold () {
	local chr="$1" outdir="$2" base
	base="${study_name}_chr${chr}_15_scaffold"
	ld_pack_archive "${section_15_upload_dir}" "${base}" "${outdir}" \
		manifest.json variants.tsv.gz D.npy B.npy checksums.json qc_report.txt
	artifact_pair_exists "${base}"
}

pack_chunk () {
	local chr="$1" outdir="$2" chunk="$3" base chunk_dir
	base="${study_name}_chr${chr}_15_chr${chr}_${chunk}"
	chunk_dir="${outdir}/A_blocks/chr${chr}/${chunk}"
	if artifact_pair_exists "${base}"; then
		rm -rf "${chunk_dir}"
		return 0
	fi
	if [ ! -d "${chunk_dir}" ]; then
		echo "Problem: raw chunk ${chunk_dir} is missing and ${base} has not been packaged" >&2
		return 1
	fi
	ld_pack_archive "${section_15_upload_dir}" "${base}" "${outdir}/A_blocks" "chr${chr}/${chunk}"
	artifact_pair_exists "${base}"
	rm -rf "${chunk_dir}"
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
	check_required_files "${outdir}" "${chr}"
	if [ -f "${outdir}/.packaged" ] && all_artifacts_exist "${chr}" "${outdir}"; then
		echo "[15b] chr${chr} already packaged; skip"
		echo "Successfully packaged LD cohort chromosome chr${chr}"
		return 0
	fi

	tmp_chunks="$(mktemp)"
	manifest_chunks "${outdir}" "${chr}" > "${tmp_chunks}"
	check_unexpected_chunks "${chr}" "${outdir}" "${tmp_chunks}"
	if [ ! -f "${outdir}/.checksums_verified" ]; then
		verify_checksums "${outdir}"
		touch "${outdir}/.checksums_verified"
	fi
	while IFS= read -r chunk; do
		pack_chunk "${chr}" "${outdir}" "${chunk}"
	done < "${tmp_chunks}"
	rm -f "${tmp_chunks}"

	pack_scaffold "${chr}" "${outdir}"
	all_artifacts_exist "${chr}" "${outdir}"
	rm -rf "${outdir}/A_blocks"
	touch "${outdir}/.packaged"
	echo "Successfully packaged LD cohort chromosome chr${chr}"
}

if [ "${ld_chromosome}" = "all" ]; then
	echo "Problem: section 15 must be packaged one chromosome at a time, immediately after 15a."
	echo "Run, for example: ld_chromosome=22 bash 15b-ld_check_compress_data.sh -c config"
	exit 1
fi

chromosome_dir="$(resolve_chromosome_dir "${ld_chromosome}")"
chromosome_log="${section_15_dir}/logs_b/chr${ld_chromosome}.log"
mkdir -p "$(dirname "${chromosome_log}")"
exec &> >(tee "${chromosome_log}")

process_chromosome "${ld_chromosome}" "${chromosome_dir}"
