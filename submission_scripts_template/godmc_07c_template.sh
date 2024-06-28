#!/bin/bash -l
#SBATCH --job-name=GoDMC_07c
#SBATCH --output=../job_reports/GoDMC_07c_%A_%a
#SBATCH --partition gpu,cpu
#SBATCH --mem=16GB
#SBATCH --ntasks=1
#SBATCH --time=6:0:0

sbatch ../07c-check_vmeQTL_results.sh drm
sbatch ../07c-check_vmeQTL_results.sh svlm
sbatch ../07c-check_vmeQTL_results.sh BF
