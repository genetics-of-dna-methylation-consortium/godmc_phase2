#!/bin/bash

./resources/setup.sh "$@"
exec &> >(tee ${section_04a_logfile})
print_version
#Please read reasources/bin/hase/README_2.md
#An expample is also provided below

mkdir -p ${hase_dir_in}
mkdir -p ${hase_converting}
cp ${bfile}.bim	${hase_dir_in}
cp ${bfile}.fam ${hase_dir_in}
cp ${bfile}.bed ${hase_dir_in}

python ${hase}/hase.py \
    -mode converting \
    -g ${hase_dir_in} \
    -o ${hase_converting} \
    -study_name ${study_name} # the name for your study 

echo "Successfully converted genetic data"
