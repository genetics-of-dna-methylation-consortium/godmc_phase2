# Hase for Meta Analysis
This is a python 3.x compatible version of Hase which was initially developed in python 2.7
To use the package, please read the wiki and follow the instructions given in the original version of the application here - 
https://github.com/roshchupkin/hase/wiki


## Setup

```bash
#Create the Conda Environment
conda env create -f environment.yml
conda activate hase
```

## Task 1: Converting
This converts the data into "genotype", "individuals" and "probes". The steps are as follows -

1. Create and specify the output folder (mkdir and -o)
2. Specify the input (the location of the plink folder) (-g)
3. Specify the name of the project (-study_name)
An example is provided below

```bash
mkdir -p ./test/converted_files

python hase.py \
-mode converting \
-g ./test/plink \
-o ./test/converted_files \
-study_name go_dmc_1 \
```

## Task 2: Mapping
This maps the converted data to the reference panel. However, the sub-task 2.1 might be necessary

### Sub-task 2.1 (Optional): Inverting the Probes
In case the allele values are flipped and there is a need to invert their positions, you can use this script provided
here with the following steps-

1. Specify the location of the probes. It will be in the probes folder of the "converting" step
output (-f)
2. Specify the studyname (-n)
3. specify the reference panel directory

```bash
python ./added/invert_probes.py \
-f ./test/converted_files/probes \
-n go_dmc_1 \
```

### Sub-task 2.2: The Mapping
1. Create and specify the output folder (mkdir and -o)
2. Specify the location of the converted files (the output of task 1) (-g)
3. Specify the name of the project (-study_name)
4. Specify the reference panel (-ref_name)
An example is provided below

```bash
mkdir -p ./test/mapped_files

python ./tools/mapper.py \
-g ./test/converted_files \
-o ./test/mapped_files \
-study_name go_dmc_1 \
-ref_name ref-hrc
```

## Task 3: Encoding
This performs encoding. The steps are as follows

1. Create and specify the output folder (mkdir and -o)
2. Specify the name of the project (-study_name)
3. Specify the location of the converted files (the output of task 1) (-g)
4. Specify the location of the mapped files (-mapper)
5. Specify the location of the phenotype data. Note that this should be a .csv file (-ph) which CpGs in rows and individuals in cols. The file should be named "methylation_data.csv"
6. Specify the reference panel (-ref_name)

An example is provided below

```bash
mkdir -p ./test/encoded_files

python hase.py \
-mode encoding \
-study_name go_dmc_1 \
-g ./test/converted_files \
-o ./test/encoded_files \
-mapper ./test/mapped_files \
-ph ./test/hase_pheno \
-ref_name ref-hrc
```

## Task 4: Single Site Analysis
This performs single site analysis. The steps are as follows

1. Create and specify the output folder (mkdir and -o)
2. Specify the name of the project (-study_name)
3. Specify the location of the converted files (the output of task 1) (-g)
4. Specify the location of the phenotype data. Note that this should be a .csv file (-ph)
5. Specify the location of the covariates data (-cov)
6. Specify the location of the mapped files (-mapper)
7. Specify the reference panel (-ref_name)

An example is provided below

```bash
mkdir -p ./test/single_site_files

python hase.py \
-mode single-meta \
-study_name go_dmc_1 \
-g ./test/converted_files \
-ph ./test/hase_pheno \
-cov ./test/hase_cov  \
-mapper ./test/mapped_files  \
-o ./test/single_site_files \
-ref_name ref-hrc
```



## Task 5: Meta Analysis
This performs single site analysis. The steps are as follows

1. Create and specify the output folder (mkdir and -o)
2. Create a separate folder for each study. In each, create a 4 folders - genotype, individuals, phenotype and probes.
3. Copy encoded results from genotype, individuals and phenotypes into the created folders in point 2 above.
4. Copy probes from results of Task 1 into probes ceated in point 2 above.
5. Create a separate folder for each single site cohort analysis and copy the results from the single site analysis in task 4 into it
6. Create a folder and copy the mapping results from task 2 of each cohort into it
7. Specify the location of the genotype (-g)
8. Specify the names of the cohorts (-study_name)
9. Specify the location of the phenotype (-ph)
10. Specify the location of the combined mapped files (-mapper)
11. specify the mode at meta-stage (-mode)
12. Specify the node, if desired (-node)
13. Speicy threshold and MAF (-th and -maf)

An example is provided below

```bash
mkdir -p ./test/meta_analysis

python hase.py \
-g ./test/meta_inputs_x/folders  ./test/meta_inputs_y/folders \
-study_name study_1 study_2 \
-ph ./test/meta_inputs_x/folders/phenotype ./test/meta_inputs_y/folders/phenotype \
-derivatives ./meta_inputs_x/derivatives ./test/meta_inputs_y/derivatives \
-mapper ./test/combined_mapped_files \
-o ./test/meta_analysis \
-mode meta-stage \
-node 100 50 \
-ref_name ref-hrc \
-th 0 \
-maf 0 \
-encoded 1 1
```