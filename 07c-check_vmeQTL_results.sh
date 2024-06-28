#!/bin/bash -l

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
exec &> >(tee ${section_07c_logfile}_$1)
print_version

# check output of OSCA

method=$1

res_count=`ls ${section_07_dir}/vQTL_${method}_cis_genetic*.besd | wc -l`
tabcount=`wc -l ${section_07_dir}/tabfile.info1|cut -d ' ' -f1`

if [ $res_count == $tabcount ]
then
   echo "using $method method to detect vmeQTLs is successful"
else
   for g_chunk in $(seq 1 ${genetic_chunks})
   do
       chr_temp=`awk 'BEGIN{FS=" "}{if($1=="'"${g_chunk}"'") print $2}' ${section_07_dir}/tabfile.info1`
       for chr in ${chr_temp}
       do
                if ! [ -f ${section_07_dir}/vQTL_${method}_cis_genetic${g_chunk}_cpgchr${chr}.besd ]
                then    
                    echo "resubmit vmeQTL detection job - vmeQTL method: ${method}, genetic chunk: ${g_chunk}, chr: ${chr}"
                    sbatch 07b-run_cis_vmeQTL.sh ${method} ${chr} ${g_chunk}
                fi
       done
   done
fi
