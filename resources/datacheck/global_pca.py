# import onnx
# import random
# from gnomad.sample_qc.ancestry import (
#     apply_onnx_classification_model,
#     assign_population_pcs,
#     pc_project,
# )
# from gnomad.utils.filtering import filter_to_adj
import hail as hl
import psutil
import pandas as pd
import matplotlib.pyplot as plt
import gzip
import json
import warnings
import sys
import os

logfile = sys.argv[1]
bfile = sys.argv[2]
study_name = sys.argv[3]
home_directory = sys.argv[4]
scripts_directory = sys.argv[5]

print("Study name:", study_name)
mem = psutil.virtual_memory()
avail_gb = mem.available / (1024**3)
driver_mem_gb = int(avail_gb * 0.8)

print(f"Available mem: {avail_gb:.2f} GB")
print(f"Total mem: {mem.total / (1024**3):.2f} GB")

genome_build = 37

slurm_cores = os.environ.get('SLURM_CPUS_ON_NODE')
if slurm_cores is not None and slurm_cores.isdigit():
    cores_to_use = int(slurm_cores)
else:
    cores_to_use = os.cpu_count()

print(f"Using {cores_to_use} cores for Hail")
hl.init(
    backend="spark",
    local=f"local[{cores_to_use}]",
    log=logfile,
    spark_conf={
        'spark.driver.memory': f'{driver_mem_gb}g',
    }
)

# import genotype data, if genome build is 37, lift over to 38
if genome_build == 37:
    print("reading genotype data")
    dat = hl.import_plink(bed=f'{bfile}.bed',
                          bim=f'{bfile}.bim', 
                          fam=f'{bfile}.fam',
                          reference_genome='GRCh37')
    chain_file = f"{scripts_directory}/resources/genetics/references_grch37_to_grch38.over.chain.gz"
    rg37 = hl.get_reference('GRCh37')
    rg38 = hl.get_reference('GRCh38')
    rg37.add_liftover(chain_file, rg38)
    dat = dat.annotate_rows(new_locus = hl.liftover(dat.locus, 'GRCh38'))
    dat = dat.filter_rows(hl.is_defined(dat.locus))
    dat = dat.key_rows_by(
        locus = dat.new_locus,
        alleles = dat.alleles)
    dat = dat.drop('new_locus')
    print("Successfully lifted over to GRCh38")
elif genome_build == 38:
    print("reading genotype data")
    dat = hl.import_plink(bed=f'{bfile}.bed',
                          bim=f'{bfile}.bim', 
                          fam=f'{bfile}.fam',
                          reference_genome='GRCh38')
    print("Original genome build: GRCh38")

# ref panel loadings
ref_pca_loadings = hl.import_table(
    f'{scripts_directory}/resources/genetics/hgdp_tgp_unrel_pass_filtGBMI_strictpruned_loadings.tsv',
    delimiter='\t',
    types={
        'locus': hl.tlocus('GRCh38'),
        'alleles': hl.tarray(hl.tstr),
        'loadings': hl.tarray(hl.tfloat64),
        'pca_af': hl.tfloat64})
ref_pca_loadings = ref_pca_loadings.key_by('locus')

# filter to locus rows that are in the ref panel
dat_filter = dat.filter_rows(hl.is_defined(ref_pca_loadings[dat.locus]))

# identify the alleles need flip to match the ref loadings
dat_filter = dat_filter.annotate_rows(
    ref_alleles = ref_pca_loadings[dat_filter.locus].alleles,
    need_flip = (dat_filter.alleles[0] == ref_pca_loadings[dat_filter.locus].alleles[1]) & 
                (dat_filter.alleles[1] == ref_pca_loadings[dat_filter.locus].alleles[0])
)

dat_filter = dat_filter.key_rows_by('locus')

dat_filter = dat_filter.annotate_rows(
    alleles = hl.if_else(
        dat_filter.need_flip,
        [dat_filter.alleles[1], dat_filter.alleles[0]],
        dat_filter.alleles
    )
)

dat_filter = dat_filter.transmute_entries(
    GT = hl.if_else(
        dat_filter.need_flip,
        hl.if_else(
            hl.is_defined(dat_filter.GT),
            hl.call(2 - dat_filter.GT[0], 2 - dat_filter.GT[1]),
            hl.missing(hl.tcall)
        ),
        dat_filter.GT
    )
)
dat_filter = dat_filter.drop('ref_alleles', 'need_flip')

dat_filter = dat_filter.key_rows_by('locus', 'alleles')
ref_pca_loadings = ref_pca_loadings.key_by('locus', 'alleles')

# find the missing rate of SNPs in loadings
nrows, ncols = dat_filter.count()
print(f"Filtered sample: {nrows} rows, {ncols} samples")

call_rate = nrows / ref_pca_loadings.count()

if call_rate >= 0.9:
    print(f"Call rate: {call_rate:.2%}")
else:
    print(f"Call rate: {call_rate:.2%}, which is less than 90%. Please contact with Haotian Tang (haotian.tang@bristol.ac.uk)")

ht_projections = hl.experimental.pc_project(dat_filter.GT, ref_pca_loadings.loadings, ref_pca_loadings.pca_af)
ht_projections = ht_projections.transmute(**{f"PC{i}": ht_projections.scores[i - 1] for i in range(1, 21)})

ht_projections = ht_projections.key_by("s")
ht_projections = ht_projections.select(
    **{"IID": hl.str(ht_projections.key)}, 
    **{f"PC{i}": ht_projections[f"PC{i}"] for i in range(1, 21)}
)

# export projections
study_projection = ht_projections.to_pandas()

# extract the first two PCs from previously calculated PCA scores
ref_score = pd.read_csv(f'{scripts_directory}/resources/genetics/hgdp_tgp_unrel_pass_filtGBMI_strictpruned_scores.tsv', sep='\t')
scores_list = ref_score['scores'].apply(lambda x: eval(x) if isinstance(x, str) else x)
ref_score['PC1'] = scores_list.apply(lambda x: x[0])
ref_score['PC2'] = scores_list.apply(lambda x: x[1])

# extract genetic region from gnomad sample summary data
with gzip.open(f'{scripts_directory}/resources/genetics/release_3.1.2_vcf_genomes_gnomad.genomes.v3.1.2.hgdp_1kg_subset_sample_meta.tsv.bgz', 'rt') as f:
    meta_df = pd.read_csv(f, sep='\t')

meta_df_filt = meta_df[meta_df['s'].isin(ref_score['s'])]
meta_subset = meta_df_filt[['s', 'hgdp_tgp_meta']]

def extract_genetic_region(hgdp_meta):
    """from hgdp_tgp_meta dictionary, extract genetic_region"""
    if pd.isna(hgdp_meta):
        return None
    try:
        if isinstance(hgdp_meta, str):
            meta_dict = json.loads(hgdp_meta)
        else:
            meta_dict = hgdp_meta
        return meta_dict.get('genetic_region', None)
    except (json.JSONDecodeError, AttributeError, TypeError):
        return None

meta_subset['genetic_region'] = meta_subset['hgdp_tgp_meta'].apply(extract_genetic_region)
print(meta_subset['genetic_region'].value_counts())

# merge the genetic region with the ref score
hgdp_tgp_ref = pd.merge(ref_score, meta_subset[['s', 'genetic_region']], on='s', how='left')

eur_subset = hgdp_tgp_ref[hgdp_tgp_ref['genetic_region'] == 'EUR']

if not eur_subset.empty:
    if eur_subset['PC1'].mean() > 0:
        print("EUR PC1 mean > 0, flipping PC1 for all samples and ref")
        hgdp_tgp_ref['PC1'] = hgdp_tgp_ref['PC1'] * -1
        study_projection['PC1'] = study_projection['PC1'] * -1
    if eur_subset['PC2'].mean() > 0:
        print("EUR PC2 mean > 0, flipping PC2 for all samples and ref")
        hgdp_tgp_ref['PC2'] = hgdp_tgp_ref['PC2'] * -1
        study_projection['PC2'] = study_projection['PC2'] * -1

# plot
plt.figure(figsize=(10, 6), constrained_layout=True)

for pop in hgdp_tgp_ref['genetic_region'].unique():
    if pd.notna(pop):
        sub = hgdp_tgp_ref[hgdp_tgp_ref['genetic_region'] == pop]
        plt.scatter(sub['PC1'], sub['PC2'], label=pop, s=10, alpha=0.5)

plt.scatter(study_projection['PC1'], study_projection['PC2'],
            color='black', marker='x', s=10, label=f'{study_name}')

plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title('Projection of Samples onto gnomAD PCA space')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

if 'raw' in logfile:
    print("Logfile indicates raw data")
    suffix = '_rawdat'
elif 'clean' in logfile:
    print("Logfile indicates clean data")
    suffix = '_cleandat'
else:
    suffix = ''

plt.savefig(f'{home_directory}/results/01/{study_name}_globalPCA{suffix}.png', dpi=300, bbox_inches='tight')