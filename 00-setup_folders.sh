#!/bin/bash

chmod +x ./resources/setup.sh
./resources/setup.sh "$@"

print_version

mkdir -p ${home_directory}/processed_data/ids
mkdir -p ${home_directory}/processed_data/genetic_data
mkdir -p ${home_directory}/processed_data/covariate_data
mkdir -p ${home_directory}/processed_data/methylation_data
mkdir -p ${home_directory}/processed_data/cellcounts
mkdir -p ${home_directory}/processed_data/inversions
mkdir -p ${home_directory}/processed_data/genetic_data/tabfile
mkdir -p ${home_directory}/job_reports
mkdir -p ${hase}/data

for i in {1..9}
do
   mkdir -p ${home_directory}/results/0$i
done

for i in {10..14}
do
   mkdir -p ${home_directory}/results/$i
done

mkdir -p "${section_01_dir}/logs/"
mkdir -p "${section_02_dir}/logs_a/"
mkdir -p "${section_02_dir}/logs_b/"
mkdir -p "${section_03_dir}/logs_a/"
mkdir -p "${section_03_dir}/logs_b/"
mkdir -p "${section_03_dir}/logs_c/"
mkdir -p "${section_03_dir}/logs_d/"
mkdir -p "${section_03_dir}/logs_e/"
mkdir -p "${section_03_dir}/logs_f/"
mkdir -p "${section_03_dir}/logs_g/"
mkdir -p "${section_04_dir}/logs_a/"
mkdir -p "${section_04_dir}/logs_b/"
mkdir -p "${section_04_dir}/logs_c/"
mkdir -p "${section_04_dir}/logs_d/"
mkdir -p "${section_04_dir}/logs_e/"
mkdir -p "${section_04_dir}/logs_f/"
mkdir -p "${section_05_dir}/logs_a/"
mkdir -p "${section_05_dir}/logs_b/"
mkdir -p "${section_05_dir}/logs_c/"
mkdir -p "${section_06_dir}/logs/"
mkdir -p "${section_07_dir}/logs_a/"
mkdir -p "${section_07_dir}/logs_b/"
mkdir -p "${section_07_dir}/logs_c/"
mkdir -p "${section_07_dir}/logs_d/"
mkdir -p "${section_08_dir}/logs_a/"
mkdir -p "${section_08_dir}/logs_b/"
mkdir -p "${section_09_dir}/logs/"
mkdir -p "${section_10_dir}/logs_a/"
mkdir -p "${section_10_dir}/logs_b/"
mkdir -p "${section_11_dir}/logs_a/"
mkdir -p "${section_11_dir}/logs_b/"
mkdir -p "${section_11_dir}/logs_c/"
mkdir -p "${section_12_dir}/logs/"
mkdir -p "${section_13_dir}/logs_a/"
mkdir -p "${section_13_dir}/logs/"
mkdir -p "${section_14_dir}/logs/"

chmod +x *.sh
chmod +x ./resources/bin/*
