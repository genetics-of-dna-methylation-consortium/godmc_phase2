#!/bin/bash

cd ..
source ./config

for i in $(seq 1 ${meth_chunks});
do 
    sbatch --mem 256G ${scripts_directory}/06f.interaction_trans_candidateSNPs.sh $i
done
