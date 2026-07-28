import pandas as pd
import torch
import tensorqtl
import sys
import pyarrow
from tensorqtl import genotypeio, cis, trans

plink_prefix_path = str(sys.argv[1])
DNAm_bed = str(sys.argv[2])
environment_file = str(sys.argv[3])
output_prefix_file = str(sys.argv[4])
maf_thres = float(sys.argv[5])
pval_thres = float(sys.argv[6])

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(DNAm_bed)

E_df = pd.read_csv(environment_file, sep='\t', index_col=0)
E_df.index = E_df.index.astype(str)
phenotype_df = phenotype_df[E_df.index]

pr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pr.load_genotypes()
variant_df = pr.bim.set_index('snp')[['chrom', 'pos']]

def runGE(Env):
    E_df_tmp = E_df[[Env]].dropna()
    interaction_series = E_df_tmp.squeeze()
    trans_df = trans.map_trans(genotype_df, phenotype_df, covariates_df=None, 
                               interaction_s=interaction_series, 
                               return_sparse=True, 
                               pval_threshold=pval_thres, maf_threshold=maf_thres)
    trans_df = trans.filter_cis(trans_df, phenotype_pos_df, variant_df, window=2000000)
    output_file = output_prefix_file+'.'+Env+'.parquet'
    trans_df.to_parquet(output_file, engine='pyarrow', compression='snappy')

for column_name in E_df.columns:
    print(f"Starting execution for: {column_name}")
    runGE(column_name)
