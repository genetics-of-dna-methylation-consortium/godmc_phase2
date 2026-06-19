#!/usr/bin/env bash
# ld_decrypt_cohort.sh — decrypt + reassemble section-15 (LD) cohort archives
# into the cohort_stats/ tree that `ld_aggregate accumulate` consumes. No network
# I/O. Mirror of ld_encrypt_cohort.sh.
# See docs/superpowers/specs/2026-06-19-section15-ld-central-decryption-design.md
#
# Usage: ld_decrypt_cohort.sh <input_dir> <output_dir> <study_name>
#   input_dir   dir holding <study>_15_*.tgz.aes + .md5sum (left untouched)
#   output_dir  where the cohort_stats/ tree is rebuilt (== accumulate --cohort-dir)
#   study_name  selects this cohort's archives by filename prefix
#
# Override the gpg binary/wrapper for testing via the GPG env var.
set -euo pipefail

if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input_dir> <output_dir> <study_name>" >&2
  exit 2
fi

# Absolutise paths: md5sum -c runs from the staging dir, so a relative input_dir
# would otherwise break. output_dir is created first so it can be resolved too.
input_dir="$(cd "$1" && pwd)"
mkdir -p "$2"
output_dir="$(cd "$2" && pwd)"
study_name="$3"
GPG="${GPG:-gpg}"

staging="${output_dir}/.staging"
mkdir -p "${staging}"

# Decrypt ${input_dir}/$1.tgz.aes -> ${staging}/$1.tgz and verify its .md5sum.
# The input .aes is never modified or removed. The caller derives the tar path as
# "${staging}/$1.tgz" itself (do NOT capture this function's stdout — md5sum -c
# prints an "OK" line there).
decrypt_and_verify () {
  local base="$1"
  local tgz="${staging}/${base}.tgz"
  "${GPG}" --output "${tgz}" --decrypt "${input_dir}/${base}.tgz.aes"
  ( cd "${staging}" && md5sum -c "${input_dir}/${base}.md5sum" )
}

# --- 1. Scaffold (small files: manifest, variants, D, B, checksums.json) ---
scaffold="${study_name}_15_scaffold"
if [ -f "${output_dir}/manifest.json" ]; then
  echo "[ld_decrypt] skip ${scaffold} (manifest.json already present)"
else
  if [ ! -f "${input_dir}/${scaffold}.tgz.aes" ]; then
    # Try to detect what study is actually present (helps diagnose study_name mismatch).
    found=""
    shopt -s nullglob
    for f in "${input_dir}"/*_15_scaffold.tgz.aes; do
      found="${found} $(basename "${f}" _15_scaffold.tgz.aes)"
    done
    shopt -u nullglob
    if [ -n "${found}" ]; then
      echo "[ld_decrypt] ERROR: scaffold archive not found for study_name '${study_name}'; found study_name(s):${found}" >&2
    else
      echo "[ld_decrypt] ERROR: scaffold archive not found: ${input_dir}/${scaffold}.tgz.aes" >&2
    fi
    exit 1
  fi
  decrypt_and_verify "${scaffold}"
  tgz="${staging}/${scaffold}.tgz"
  tmp="${staging}/scaffold_extract.$$"
  rm -rf "${tmp}"; mkdir -p "${tmp}"
  tar xzf "${tgz}" -C "${tmp}"
  # Move every scaffold file into place, manifest.json LAST so its presence is a
  # reliable "scaffold complete" resume marker.
  for f in "${tmp}"/*; do
    [ "$(basename "${f}")" = "manifest.json" ] && continue
    mv -f "${f}" "${output_dir}/"
  done
  mv -f "${tmp}/manifest.json" "${output_dir}/manifest.json"
  rm -rf "${tmp}" "${tgz}"
  echo "[ld_decrypt] restored ${scaffold}"
fi

# --- 2. Study-name cross-check (hard fail on mismatch) ---
got_study="$(python -c 'import json,sys; print(json.load(open(sys.argv[1])).get("study_name",""))' \
  "${output_dir}/manifest.json")"
if [ "${got_study}" != "${study_name}" ]; then
  echo "[ld_decrypt] ERROR: manifest study_name '${got_study}' != requested '${study_name}'" >&2
  exit 1
fi

# --- 3. Per-chunk A_blocks archives ---
mkdir -p "${output_dir}/A_blocks"
n_chunks=0
shopt -s nullglob
for aes in "${input_dir}/${study_name}_15_chr"*"_chunk_"*".tgz.aes"; do
  base="$(basename "${aes}" .tgz.aes)"           # e.g. study_15_chr1_chunk_0
  suffix="${base#${study_name}_15_}"             # chr1_chunk_0
  chr_name="${suffix%%_chunk_*}"                 # chr1
  chunk_name="chunk_${suffix#*_chunk_}"          # chunk_0
  target="${output_dir}/A_blocks/${chr_name}/${chunk_name}"
  n_chunks=$((n_chunks + 1))
  if [ -d "${target}" ]; then
    echo "[ld_decrypt] skip ${base} (already restored)"
    continue
  fi
  decrypt_and_verify "${base}"
  tgz="${staging}/${base}.tgz"
  tmp="${staging}/chunk_extract.$$"
  rm -rf "${tmp}"; mkdir -p "${tmp}"
  tar xzf "${tgz}" -C "${tmp}"                    # tmp now holds chr<C>/chunk_<N>/...
  mkdir -p "${output_dir}/A_blocks/${chr_name}"
  mv "${tmp}/${chr_name}/${chunk_name}" "${target}"   # atomic publish
  rm -rf "${tmp}" "${tgz}"
  echo "[ld_decrypt] restored ${base}"
done
shopt -u nullglob

if [ "${n_chunks}" -eq 0 ]; then
  echo "[ld_decrypt] ERROR: no A_blocks chunk archives found for '${study_name}' in ${input_dir}" >&2
  exit 1
fi

# Best-effort integrity guard: compare against the manifest's declared chunk count.
manifest="${output_dir}/manifest.json"
if command -v jq >/dev/null 2>&1 && [ -f "${manifest}" ]; then
  expected="$(jq '[.A_blocks.chromosomes[].n_chunks] | add // 0' "${manifest}" 2>/dev/null || true)"
  if [ -n "${expected}" ] && [ "${expected}" -gt 0 ] && [ "${expected}" != "${n_chunks}" ]; then
    echo "[ld_decrypt] WARNING: manifest expects ${expected} chunks but restored ${n_chunks}" >&2
  fi
fi

rmdir "${staging}" 2>/dev/null || true
echo "[ld_decrypt] done: ${n_chunks} chunk archive(s) + scaffold in ${output_dir}"
