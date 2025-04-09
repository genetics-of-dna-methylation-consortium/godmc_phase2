#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04a_logfile})
print_version
#Please read resources/bin/hase/README_2.md
#An expample is also provided below

mkdir -p ${hase_dir_in}
mkdir -p ${hase_converting}

nX=`grep ^X ${bfile}.bim | wc -l`
if [ "$nX" -gt "0" ]
then

#perl -pe 's/^X\tX/23\t23/g' < ${bfile}.bim >${hase_dir_in}/data.bim
${plink2} --bfile ${bfile} --make-bed --output-chr 26 --out ${hase_dir_in}/data
cp ${bfile}.fam ${hase_dir_in}
cp ${bfile}.bed ${hase_dir_in}
rm ${hase_dir_in}/data.log
else
cp ${bfile}.bim ${hase_dir_in}
cp ${bfile}.fam ${hase_dir_in}
cp ${bfile}.bed ${hase_dir_in}
fi

nX=`grep ^X ${bfile}.bim | wc -l`


python ${hase}/hase.py \
    -mode converting \
    -g ${hase_dir_in} \
    -o ${hase_converting} \
    -study_name ${study_name} # the name for your study 

echo "Successfully converted genetic data"
