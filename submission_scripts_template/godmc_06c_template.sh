#!/bin/bash

## observe the chromosome information of each tab file

cd ..
source config

for i in $(seq 1 ${genetic_chunks})
do
    chr=`cut -f 1 ${tabfile}.tab.$i.bim | sort | uniq`
    echo $i $chr >> ${section_06_dir}/tabfile.GEI.info
done

count=`awk 'BEGIN{FS=" "}{print NF}' ${section_06_dir}/tabfile.GEI.info | sort | uniq`
for i in $count
do
    cut -d ' ' -f 1,$i ${section_06_dir}/tabfile.GEI.info >> ${section_06_dir}/tabfile.GEI.temp
done

awk 'BEGIN{FS=" "}{if($2>0 && $2<23) print}' ${section_06_dir}/tabfile.GEI.temp | sort -k1,1n | uniq > ${section_06_dir}/tabfile.GEI.info1

rm ${section_06_dir}/tabfile.GEI.info
rm ${section_06_dir}/tabfile.GEI.temp

cat ${section_06_dir}/tabfile.GEI.info1 | while read line; 
do 
    genetic_chunk=`echo $line | cut -d " " -f 1`
    chr=`echo $line | cut -d " " -f 2`
    sbatch --mem 64G 06c-interaction_cis.sh $genetic_chunk $chr
done
