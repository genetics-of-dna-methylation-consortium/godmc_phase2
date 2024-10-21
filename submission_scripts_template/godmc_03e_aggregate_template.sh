#!/bin/bash -l
#SBATCH --job-name=GoDMC_03e_aggregate
#SBATCH --output=../job_reports/GoDMC_03e_aggregate_%j
#SBATCH --partition gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0

cd ..
source resources/setup.sh "$@"
set -- $concatenated

bash ./resources/methylation/aggregate_adjustment2.sh
