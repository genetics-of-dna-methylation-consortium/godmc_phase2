import pandas as pd
import numpy as np
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
print("Any -9 before imputation?", (genotype_df == -9).values.any())
genotype_df = genotype_df.astype(np.float32)
genotype_df[genotype_df == -9] = np.nan
genotype_df = genotype_df.apply(lambda row: row.fillna(row.mean()), axis=1)
print("Any -9 after imputation?", (genotype_df == -9).values.any())
print("Any remaining NaNs after imputation?", genotype_df.isna().values.any())

variant_df = pr.bim.set_index('snp')[['chrom', 'pos', 'a0', 'a1']]
variant_df.to_csv(output_prefix_file + '.Allele_info.csv', index=True)

def runGE(Env, genotype_df, phenotype_df, E_df):
    E_df_tmp = E_df[[Env]].dropna()
    interaction_series = E_df_tmp.squeeze()
    if E_df_tmp[Env].std() == 0 or np.isclose(E_df_tmp[Env].std(), 0):
        print(f"[SKIP] '{Env}': Zero variance (all values are identical)")
        return

    if len(E_df_tmp) < 50:
        print(f"[SKIP] '{Env}': Too few valid samples remaining ({len(E_df_tmp)})")
        return

    E_df_tmp.index = E_df_tmp.index.astype(str)
    genotype_df.columns = genotype_df.columns.astype(str)
    common_iids = E_df_tmp.index.intersection(genotype_df.columns)

    E_df_tmp1 = E_df_tmp.loc[common_iids]
    interaction_series = E_df_tmp1.squeeze()
    genotype_df1 = genotype_df[common_iids]
    phenotype_df1 = phenotype_df[common_iids]

    same_iids_and_order = (E_df_tmp1.index.equals(genotype_df1.columns) and genotype_df1.columns.equals(phenotype_df1.columns))

    if same_iids_and_order==False:
        print("Different sample sizes in Env, Genotype, and Phenotype dataset")
        return

    try:
        trans_df = trans.map_trans(genotype_df1, phenotype_df1, covariates_df=None, 
                                   interaction_s=interaction_series, 
                                   return_sparse=True, 
                                   pval_threshold=pval_thres, maf_threshold=maf_thres)
        trans_df = trans.filter_cis(trans_df, phenotype_pos_df, variant_df, window=2000000)
        output_file = output_prefix_file+'.'+Env+'.parquet'
        trans_df.to_parquet(output_file, engine='pyarrow', compression='snappy')
        print(f"[SUCCESS] Completed GxE for: {Env}")
    except Exception as err:
        print(f"[ERROR] Failed GxE for '{Env}'. Reason: {err}")

for column_name in E_df.columns:
    print(f"Starting execution for: {column_name}")
    runGE(column_name, genotype_df, phenotype_df, E_df)
