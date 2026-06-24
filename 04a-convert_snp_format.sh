#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04a_logfile})
print_version
#Please read resources/bin/hase/README_2.md
#An expample is also provided below

mkdir -p ${light_hase_dir_in}
mkdir -p ${light_hase_converting}
mkdir -p ${light_hase}/data

# assuming ref is in ori hase folder
cp ${hase}/data/ref-hrc.ref.gz ${light_hase}/data/ref-hrc.ref.gz
cp ${hase}/data/ref-hrc.ref_info.h5 ${light_hase}/data/ref-hrc.ref_info.h5

echo "Flipping alleles into hrc reference allele order"

hrc_ref_allele="${light_hase}/data/hrc_ref_allele.txt"

zcat ${light_hase}/data/ref-hrc.ref.gz \
    | awk 'NR>1 {print $1 "\t" $3}' \
    > ${hrc_ref_allele}

echo "ref-hrc.ref.gz lines including header:"
zcat ${light_hase}/data/ref-hrc.ref.gz | wc -l

echo "${hrc_ref_allele} lines without header:"
wc -l ${hrc_ref_allele}

${plink2} \
    --bfile "${bfile}" \
    --new-id-max-allele-len 100 \
    --sort-vars \
    --set-all-var-ids @:#_\$1_\$2 \
    --ref-allele force ${hrc_ref_allele} 2 1 \
    --make-bed \
    --output-chr 26 \
    --out "${bfile}_haseinput" \
	--threads "${nthreads}"

nX=`grep ^X ${bfile}_haseinput.bim | wc -l`

if [ "$nX" -gt "0" ]
then

#perl -pe 's/^X\tX/23\t23/g' < ${bfile}.bim >${hase_dir_in}/data.bim
${plink2} --bfile ${bfile}_haseinput --make-bed --output-chr 26 --out ${hase_dir_in}/data
rm ${light_hase_dir_in}/data.log
echo "This condition should not occur anymore; but if it does, the code above will fix the chrX coding in the bim file."
else
cp ${bfile}_haseinput.bim ${light_hase_dir_in}/data.bim
cp ${bfile}_haseinput.fam ${light_hase_dir_in}/data.fam
cp ${bfile}_haseinput.bed ${light_hase_dir_in}/data.bed
fi

nX=`grep ^X ${light_hase_dir_in}/data.bim | wc -l`
if [ "$nX" -gt "0" ]
then
echo "ERROR: wrong chrX coding"
fi


python ${light_hase}/light_hase.py \
    -mode converting \
    -g ${light_hase_dir_in} \
    -o ${light_hase_converting} \
    -study_name ${study_name} # the name for your study 

echo "Successfully converted genetic data"
