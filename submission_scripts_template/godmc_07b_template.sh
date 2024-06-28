#!/bin/bash
#SBATCH --job-name=GoDMC_07b
#SBATCH --output=../job_reports/GoDMC_07b_%j
#SBATCH --partition gpu,cpu
#SBATCH --mem=16GB
#SBATCH --ntasks=1
#SBATCH --time=6:0:0

## observe the chromosome information of each tab file

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
exec &> >(tee ${section_07b_logfile})
print_version

for i in $(seq 1 ${genetic_chunks})
do
chr=`cut -f 1 ${tabfile}.tab.$i.bim | sort | uniq`
echo $i $chr >> ${section_07_dir}/tabfile.info
done

count=`awk 'BEGIN{FS=" "}{print NF}' ${section_07_dir}/tabfile.info | sort | uniq`
for i in $count
do
cut -d ' ' -f 1,$i ${section_07_dir}/tabfile.info >> ${section_07_dir}/tabfile.temp
done

awk 'BEGIN{FS=" "}{if($2>0) print}' ${section_07_dir}/tabfile.temp | sort -k1,1n | uniq > ${section_07_dir}/tabfile.info1

cat ${section_07_dir}/tabfile.info1 | while read line; 
do 
genetic_chunk=`echo $line | cut -d " " -f 1`
chr=`echo $line | cut -d " " -f 2`
sbatch 07b-run_cis_vmeQTL.sh drm $chr $genetic_chunk
sbatch 07b-run_cis_vmeQTL.sh svlm $chr $genetic_chunk
sbatch 07b-run_cis_vmeQTL.sh BF $chr $genetic_chunk
done
