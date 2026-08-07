#!/bin/bash

cd ..
source ./config

for i in $(seq 1 22);
do 
    sbatch --mem 128G ${scripts_directory}/06g-interaction_trans_candidateCpGs.sh $i
done
