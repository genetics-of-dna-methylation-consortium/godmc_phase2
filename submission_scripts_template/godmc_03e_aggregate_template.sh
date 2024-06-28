#!/bin/bash -l
#SBATCH --job-name=GoDMC_03e_aggregate
#SBATCH --output=../job_reports/GoDMC_03e_aggregate_%j
#SBATCH --partition gpu,cpu
#SBATCH --mem=64GB
#SBATCH --ntasks=8
#SBATCH --time=6:0:0

cd ..
set -e

# Initialize variables
config_file="./config"

# Parse options using getopts
while getopts "c:" opt; do
    case $opt in
        c) config_file=$OPTARG ;;
        *) echo "Usage: $0 -c <config_file>"
           exit 1 ;;
    esac
done

# Shift option arguments, so $1 becomes the first positional argument
shift $((OPTIND - 1))

set -e
echo "-----------------------------------------------"
echo ""
echo "Using config located at:" ${config_file}
echo ""
echo "-----------------------------------------------"
	
source ${config_file}

bash ./resources/methylation/aggregate_adjustment2.sh
