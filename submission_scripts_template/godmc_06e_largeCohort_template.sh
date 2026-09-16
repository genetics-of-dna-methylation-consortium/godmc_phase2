#!/bin/bash

cd ..
source ./config

for i in $(seq 1 22);
do 
    for chunk in $(seq 1 ${prune_sub});
    do
        sbatch --mem 128G ${scripts_directory}/06e-vmeQTL_detection_candidateCpGs.sh BF $i ${chunk}
    done
done
