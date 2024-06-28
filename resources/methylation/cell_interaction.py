import pandas as pd
import torch
import tensorqtl
import rpy2
import os
import sys
from tensorqtl import genotypeio, cis, trans, post

# define variables
genotype_tab = sys.argv[1]
methylation_file = sys.argv[2]
interaction_file = sys.argv[3]
chromosome = str(sys.argv[4])
cell_type = sys.argv[5]
output_dir = sys.argv[6]

#genotype_tab = 'processed_data/genetic_data/tabfile/data'
#methylation_file = 'processed_data/methylation_data/cell_interaction'
#interaction_file = 'processed_data/cellcounts/cellcounts.covariates.txt'
#chromosome = str(1)
#output_dir = '/scratch/prj/dtr/Groups_WorkSpace/JordanaBell/epigenetics/Analysis/subprojects/Xiaopu/GoDMC_2/EPIC_mQTL/results/06/'

plink_prefix_path = genotype_tab+'.chr'+chromosome
methylation_bed = methylation_file+'_chr'+chromosome+'.bed.gz'
interaction_factor = interaction_file
prefix = 'cell_interaction'

# load phenotypes and covariates
phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(methylation_bed)
interaction_df = pd.read_csv(interaction_factor, sep=' ', index_col=0)

# PLINK reader for genotypes
pgr = genotypeio.PlinkReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.bim.set_index('snp')[['chrom', 'pos']]

interaction_df_sub = interaction_df[[cell_type]]
phenotype_df.columns = phenotype_df.columns.astype('str')
interaction_df_sub.index = interaction_df_sub.index.astype('str')
cis.map_nominal(genotype_df=genotype_df, variant_df=variant_df,
                phenotype_df=phenotype_df,
                phenotype_pos_df=phenotype_pos_df, prefix=prefix+"_"+cell_type,
                covariates_df=None,
                interaction_df=interaction_df_sub, maf_threshold_interaction=0.05,
                run_eigenmt=True, output_dir=output_dir, write_top=True, write_stats=True)

