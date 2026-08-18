#!/bin/bash

cd ..
source ./config

for i in $(seq 1 9);
do 
    for chunk in $(seq 1 ${prune_sub});
    do
        sbatch -p interruptible_cpu,cpu --mem 128G ${scripts_directory}/06g-interaction_trans_candidateCpGs.sh $i ${chunk}
    done
done
