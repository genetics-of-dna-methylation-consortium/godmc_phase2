#!/usr/bin/env bash
# ld_encrypt_cohort.sh — tar+md5+GPG a section-15 (LD) per-chromosome cohort
# output into upload artifacts. No network I/O.
#
# This is the non-streaming sibling of the per-chunk packing in
# 15c-ld_run_upload.sh: both drive resources/genetics/ld_pack.sh, so the
# tar/md5/encrypt convention lives in exactly one place, and both emit the
# per-chromosome archive names that ld_reassemble_cohort.py consumes.
# See docs/superpowers/specs/2026-06-19-section15-ld-gpg-upload-design.md
#
# Usage: ld_encrypt_cohort.sh <cohort_stats_dir> <output_dir> <study_name>
#   cohort_stats_dir  a single-chromosome 15a --output-dir (exactly one
#                     A_blocks/chr<C> subtree, as 15a/15c produce per chromosome)
#   output_dir        where <study>_chr<C>_15_*.tgz.aes + .md5sum are written
#   study_name        cohort identifier used in artifact names
#
# Emits:
#   <study>_chr<C>_15_scaffold            (manifest/variants/D/B/checksums/qc)
#   <study>_chr<C>_15_chr<C>_chunk_<N>    (one per A_blocks chunk)
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

# Shared tar+md5+GPG primitive (single source of truth for section-15 packing).
source "$(dirname "${BASH_SOURCE[0]}")/ld_pack.sh"

ablocks_dir="${cohort_stats_dir}/A_blocks"

# --- Determine the single chromosome this cohort dir holds ---
shopt -s nullglob
chrom_dirs=("${ablocks_dir}"/chr*/)
shopt -u nullglob
if [ "${#chrom_dirs[@]}" -ne 1 ]; then
  echo "[ld_encrypt] ERROR: expected exactly one chromosome directory under ${ablocks_dir}, found ${#chrom_dirs[@]}" >&2
  exit 1
fi
chr_name="$(basename "${chrom_dirs[0]}")"   # e.g. chr22

mkdir -p "${output_dir}"

# --- 1. Scaffold bundle (small files) ---
scaffold="${study_name}_${chr_name}_15_scaffold"
ld_pack_archive "${output_dir}" "${scaffold}" "${cohort_stats_dir}" \
  manifest.json variants.tsv.gz D.npy B.npy checksums.json qc_report.txt
echo "[ld_encrypt] encrypted ${scaffold}"

# --- 2. Per-chunk A_blocks archives ---
n_chunks=0
shopt -s nullglob
for chunk_dir in "${ablocks_dir}/${chr_name}"/chunk_*; do
  [ -d "${chunk_dir}" ] || continue
  n_chunks=$((n_chunks + 1))
  chunk_name="$(basename "${chunk_dir}")"                          # e.g. chunk_0
  base="${study_name}_${chr_name}_15_${chr_name}_${chunk_name}"    # study_chr22_15_chr22_chunk_0
  ld_pack_archive "${output_dir}" "${base}" "${ablocks_dir}" "${chr_name}/${chunk_name}"
  echo "[ld_encrypt] encrypted ${base}"
done
shopt -u nullglob

if [ "${n_chunks}" -eq 0 ]; then
  echo "[ld_encrypt] ERROR: no A_blocks chunks found under ${ablocks_dir}/${chr_name}" >&2
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

echo "[ld_encrypt] done: ${n_chunks} chunk archive(s) + scaffold for ${chr_name} in ${output_dir}"
