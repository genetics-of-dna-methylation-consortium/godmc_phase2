#!/usr/bin/env bash
set -euo pipefail
# Larger-chromosome scale pilot: chr1 at PRODUCTION A-block defaults
# (block_size=4096, chunk_rows=50000, max_dense_gb=1.0). Cores/memory held at the
# chr22 baseline (8 / 32 GB) so the size/runtime scaling is comparable on this host.
# Activate the conda env in a host-portable way. Override either via env var:
#   CONDA_BASE=/path/to/miniforge3 CONDA_ENV=hail_env ./run_chr1_pilot.sh
# Otherwise the conda base is discovered from whichever `conda` is on PATH.
CONDA_ENV="${CONDA_ENV:-hail_env}"
if [ -z "${CONDA_BASE:-}" ]; then
  if command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
  else
    echo "ERROR: conda not found on PATH and CONDA_BASE not set." >&2
    echo "       Install conda/miniforge or set CONDA_BASE to its install prefix." >&2
    exit 1
  fi
fi
# shellcheck source=/dev/null
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV}"

CHR=1
OUTDIR="pilot_data/1kg_chr${CHR}"
THREADS="${THREADS:-8}"
mkdir -p "${OUTDIR}"

VCF_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
PANEL_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
VCF="${OUTDIR}/1kg.chr${CHR}.vcf.gz"
PANEL="${OUTDIR}/integrated_call_samples_v3.20130502.ALL.panel"
BFILE="${OUTDIR}/data"
COVARS="${OUTDIR}/covariates_intersectids.txt"

echo "==================== [chr${CHR}] DOWNLOAD $(date '+%H:%M:%S') ===================="
[ -f "${VCF}" ]   || wget -q -O "${VCF}" "${VCF_URL}"
[ -f "${PANEL}" ] || wget -q -O "${PANEL}" "${PANEL_URL}"

echo "==================== [chr${CHR}] PLINK CONVERT $(date '+%H:%M:%S') ===================="
if [ ! -f "${BFILE}.bed" ]; then
  plink2 \
    --vcf "${VCF}" \
    --chr ${CHR} \
    --snps-only just-acgt \
    --max-alleles 2 \
    --set-all-var-ids '@:#:$r:$a' \
    --make-bed \
    --threads "${THREADS}" \
    --out "${BFILE}"
fi

echo "==================== [chr${CHR}] COVARIATES $(date '+%H:%M:%S') ===================="
python - "${BFILE}.fam" "${PANEL}" "${COVARS}" <<'PY'
import sys
from pathlib import Path
fam_path, panel_path, out_path = map(Path, sys.argv[1:4])
sex_by_sample = {}
with panel_path.open() as handle:
    header = handle.readline().strip().split()
    sample_i = header.index("sample")
    gender_i = header.index("gender")
    for line in handle:
        row = line.strip().split()
        if not row:
            continue
        g = row[gender_i].lower()
        sex_by_sample[row[sample_i]] = "M" if g == "male" else "F"
with fam_path.open() as fam, out_path.open("w") as out:
    out.write("IID Age_numeric Sex_factor\n")
    for i, line in enumerate(fam, start=1):
        iid = line.split()[1]
        out.write(f"{iid} {40 + (i % 25)} {sex_by_sample.get(iid, 'F')}\n")
PY

echo "==================== [chr${CHR}] RUN 15a (production defaults) $(date '+%H:%M:%S') ===================="
SECTION15_DIR="${OUTDIR}/section15_chr${CHR}"
mkdir -p "${SECTION15_DIR}"
/usr/bin/time -v python resources/genetics/ld_prepare_stats.py \
  --study-name 1kg_chr${CHR}_pilot \
  --bfile "${BFILE}" \
  --covariates "${COVARS}" \
  --output-dir "${SECTION15_DIR}" \
  --log-file "${OUTDIR}/section15_chr${CHR}.log" \
  --hail-local-cores 8 \
  --a-block-size 4096 \
  --a-chunk-rows 50000 \
  --a-max-dense-gb 1.0 \
  --chromosome ${CHR}

echo "==================== [chr${CHR}] DONE $(date '+%H:%M:%S') ===================="
