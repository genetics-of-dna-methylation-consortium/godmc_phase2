#!/bin/bash -l
#SBATCH --job-name=GoDMC_07d
#SBATCH --output=job_reports/GoDMC_07d_%A_%a
#SBATCH --partition interruptible_cpu
#SBATCH --mem=16GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0
#SBATCH --array=1-22

sbatch ../07d-tar_results.sh ${SLURM_ARRAY_TASK_ID}
