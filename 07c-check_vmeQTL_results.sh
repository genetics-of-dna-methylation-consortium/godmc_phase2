#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_07c_logfile}_$1)
print_version

# check output of OSCA

method=$1

for g_chunk in $(seq 1 ${genetic_chunks}); do
    chr_temp=`awk 'BEGIN{FS=" "}{if($1=="'"${g_chunk}"'") print $2}' ${section_07_dir}/tabfile.info1`
    for chr in ${chr_temp}; do

    ### remove the old results with drm and svlm name
    ### the old files look like vQTL_drm_cis_genetic${g_chunk}_cpgchr${i}.besd, while the new file format is vQTL_drm_cis_genetic${g_chunk}_cpgchr${i}_1_1.besd
    if [ "${method}" == "drm" ] || [ "${method}" == "svlm" ]; then
        if [ -f "${section_07_dir}/vQTL_${method}_cis_genetic${g_chunk}_cpgchr${chr}.besd" ]; then
            rm "${section_07_dir}/vQTL_${method}_cis_genetic${g_chunk}_cpgchr${chr}."*
        fi
    fi

    if [[ "$chr" =~ ^[0-9]+$ ]]; then
        if  [ "$chr" -le 22 ]; then
                if [ ${method} = "BF" ]; then
                    if grep -q "Analysis finished" ${section_07_dir}/vQTL_BF_cis_genetic${g_chunk}_cpgchr${chr}_1_1.log; then
                        echo "using $method method to detect vmeQTLs of genetic chunk: ${g_chunk} and chr: ${chr} is successful"
                    else
                        echo "resubmit vmeQTL detection job - vmeQTL method: ${method}, genetic chunk: ${g_chunk}, chr: ${chr}"
                        sbatch 07b-run_cis_vmeQTL.sh ${method} ${chr} ${g_chunk}
                    fi
                fi

                if [ ${method} = "drm" ]; then
                    if grep -q "besd file was writen." ${section_07_dir}/vQTL_drm_cis_genetic${g_chunk}_cpgchr${chr}_1_1.log; then
                        echo "using $method method to detect vmeQTLs of genetic chunk: ${g_chunk} and chr: ${chr} is successful"
                    else
                        echo "resubmit vmeQTL detection job - vmeQTL method: ${method}, genetic chunk: ${g_chunk}, chr: ${chr}"
                        sbatch 07b-run_cis_vmeQTL.sh ${method} ${chr} ${g_chunk}
                    fi
                fi

                if [ ${method} = "svlm" ]; then
                    if grep -q "besd file was writen." ${section_07_dir}/vQTL_svlm_cis_genetic${g_chunk}_cpgchr${chr}_1_1.log; then
                        echo "using $method method to detect vmeQTLs of genetic chunk: ${g_chunk} and chr: ${chr} is successful"
                    else
                        echo "resubmit vmeQTL detection job - vmeQTL method: ${method}, genetic chunk: ${g_chunk}, chr: ${chr}"
                        sbatch 07b-run_cis_vmeQTL.sh ${method} ${chr} ${g_chunk}
                    fi
                fi
        fi
    fi
    done
done
