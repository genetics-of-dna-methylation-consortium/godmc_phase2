#!/usr/bin/env bash
# ld_encrypt_cohort.sh — stage + symmetric-GPG-encrypt section-15 (LD) cohort
# outputs for upload. No network I/O.
# See docs/superpowers/specs/2026-06-19-section15-ld-gpg-upload-design.md
#
# Usage: ld_encrypt_cohort.sh <cohort_stats_dir> <output_dir> <study_name>
#   cohort_stats_dir  the 15a --output-dir (e.g. results/15/cohort_stats)
#   output_dir        where .tgz.aes + .md5sum are written (created if absent)
#   study_name        cohort identifier used in artifact names
#
# Override the gpg binary/wrapper for testing via the GPG env var.
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <cohort_stats_dir> <output_dir> <study_name>" >&2
  exit 2
fi

cohort_stats_dir="$1"
output_dir="$2"
study_name="$3"
GPG="${GPG:-gpg}"

# $1 = archive basename (no extension); ${output_dir}/$1.tgz must already exist.
# Writes ${1}.md5sum (of the plaintext tar) and ${1}.tgz.aes, then removes the tar.
stage_archive () {
  local base="$1"
  ( cd "${output_dir}" && md5sum "${base}.tgz" > "${base}.md5sum" )
  "${GPG}" --output "${output_dir}/${base}.tgz.aes" \
    --symmetric --cipher-algo AES256 "${output_dir}/${base}.tgz"
  rm -f "${output_dir}/${base}.tgz"
}

mkdir -p "${output_dir}"

# $1 = archive basename; returns 0 (already done) if .tgz.aes and .md5sum exist.
already_done () {
  [ -f "${output_dir}/$1.tgz.aes" ] && [ -f "${output_dir}/$1.md5sum" ]
}

# --- 1. Scaffold bundle (small files) ---
scaffold="${study_name}_15_scaffold"
if already_done "${scaffold}"; then
  echo "[ld_encrypt] skip ${scaffold} (already encrypted)"
else
  tar czf "${output_dir}/${scaffold}.tgz" -C "${cohort_stats_dir}" \
    manifest.json variants.tsv.gz D.npy B.npy checksums.json
  stage_archive "${scaffold}"
  echo "[ld_encrypt] encrypted ${scaffold}"
fi

# --- 2. Per-chunk A_blocks archives ---
ablocks_dir="${cohort_stats_dir}/A_blocks"
n_chunks=0
shopt -s nullglob
for chunk_dir in "${ablocks_dir}"/chr*/chunk_*; do
  [ -d "${chunk_dir}" ] || continue
  n_chunks=$((n_chunks + 1))
  chr_name="$(basename "$(dirname "${chunk_dir}")")"   # e.g. chr1
  chunk_name="$(basename "${chunk_dir}")"              # e.g. chunk_0
  base="${study_name}_15_${chr_name}_${chunk_name}"    # e.g. study_15_chr1_chunk_0
  if already_done "${base}"; then
    echo "[ld_encrypt] skip ${base} (already encrypted)"
    continue
  fi
  tar czf "${output_dir}/${base}.tgz" -C "${ablocks_dir}" "${chr_name}/${chunk_name}"
  stage_archive "${base}"
  echo "[ld_encrypt] encrypted ${base}"
done
shopt -u nullglob

if [ "${n_chunks}" -eq 0 ]; then
  echo "[ld_encrypt] ERROR: no A_blocks chunks found under ${ablocks_dir}" >&2
  exit 1
fi

# Best-effort integrity guard: compare against the manifest's declared chunk count.
manifest="${cohort_stats_dir}/manifest.json"
if command -v jq >/dev/null 2>&1 && [ -f "${manifest}" ]; then
  expected="$(jq '[.A_blocks.chromosomes[].n_chunks] | add // 0' "${manifest}" 2>/dev/null || true)"
  if [ -n "${expected}" ] && [ "${expected}" -gt 0 ] && [ "${expected}" != "${n_chunks}" ]; then
    echo "[ld_encrypt] WARNING: manifest expects ${expected} chunks but found ${n_chunks}" >&2
  fi
fi

echo "[ld_encrypt] done: ${n_chunks} chunk archive(s) + scaffold in ${output_dir}"
