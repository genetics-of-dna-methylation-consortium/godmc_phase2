#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

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
                if ! [ -f ${section_07_dir}/vQTL_${method}_cis_genetic${g_chunk}_cpgchr${chr}.besd ] && [[ "$chr" =~ ^[0-9]+$ ]]
                then
			if  [ "$chr" -le 22 ]
			then   
                    	echo "resubmit vmeQTL detection job - vmeQTL method: ${method}, genetic chunk: ${g_chunk}, chr: ${chr}"
	            	sbatch 07b-run_cis_vmeQTL.sh ${method} ${chr} ${g_chunk}
			fi
                fi
       done
   done
fi
