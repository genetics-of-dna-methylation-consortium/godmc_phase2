import pandas as pd
import torch
import tensorqtl
import sys
import pickle
from tensorqtl import genotypeio, cis, trans

list1 = str(sys.argv[1])
plink_prefix_path = str(sys.argv[2])
DNAm_bed = str(sys.argv[3])
chunk = str(sys.argv[4])
chrom = str(sys.argv[5])
environment_file = str(sys.argv[6])
output_dir = str(sys.argv[7])

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(DNAm_bed)
phenotype_pos_df['chr'] = phenotype_pos_df['chr'].str.replace('chr', '', regex=False)

E_df = pd.read_csv(environment_file, sep='\t', index_col=0)
E_df.index = E_df.index.astype(str)
phenotype_df = phenotype_df[E_df.index]

pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

def runGE(Env):
    E_df_tmp = E_df[[Env]].dropna()
    prefix_out = "GEI_chunk"+str(chunk)+"_chr"+str(chrom)+"_E_"+Env+".candidate"
    mapping_df = pd.read_csv(list1, sep="\t", names=['SNP','CpG','Pair'], header=None)
    cis_df = cis.map_nominal(mapping_df, genotype_df, variant_df, 
                phenotype_df, phenotype_pos_df, prefix_out,
                covariates_df=None,
                interaction_df=E_df_tmp, maf_threshold_interaction=0.001,
                run_eigenmt=False, output_dir=output_dir, write_top=False, write_stats=True)

for column_name in E_df.columns:
    print(f"Starting execution for: {column_name}")
    runGE(column_name)
