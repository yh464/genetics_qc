import pandas as pd
import os
pd.options.mode.chained_assignment = None
ref = pd.read_table('/rds/project/rds-Nl99R8pHODQ/toolbox/magma/ENSG.gene.loc', header = None, index_col = 0)
cepo_score = pd.read_table('/rds/project/rds-Nl99R8pHODQ/multiomics/gene_score/wang_2025_neocortex.cepo.txt', index_col = 'gene')
topdecile = open('/rds/project/rds-Nl99R8pHODQ/multiomics/gene_set/wang_2025_topdecile.txt').read().splitlines()
eregulons = open('/rds/project/rds-Nl99R8pHODQ/multiomics/gene_set/wang_2025_eregulons.txt').read().splitlines()
ref.index.name = None
ref.columns = ['CHR','FROM','TO','DIR','LABEL']
# flanking region 10 kb to be consistent with mixer
ref['FROM'] -= 10000; ref['TO'] += 10000
ref = ref.loc[~ref.CHR.isin(['X','Y']),:]
ref.CHR = ref.CHR.astype(int)
cepo_score = cepo_score.loc[:, cepo_score.columns.str.contains('Type')] # only care about cell type
cepo_score.columns = cepo_score.columns.str.replace('.Type','')
all_scores = cepo_score
for line in topdecile:
    line = line.split()
    all_scores.loc[:,f'topdecile.{line[0]}'] = 0
    all_scores.loc[[x for x in line[1:] if x in all_scores.index], f'topdecile.{line[0]}'] = 1
for line in eregulons:
    line = line.split()
    all_scores.loc[:,f'eregulons.{line[0]}'] = 0
    all_scores.loc[[x for x in line[1:] if x in all_scores.index], f'eregulons.{line[0]}'] = 1
all_scores = all_scores.copy()
print('Loaded scores')

import sys
from time import perf_counter as t
print(sys.argv)
chrom = int(sys.argv[1])
snps_chr = pd.read_table(f'/rds/project/rds-Nl99R8pHODQ/ref/1000g_eur_ldsc/chr{chrom}.bim', header = None, 
    usecols = [0,1,3])
snps_chr.columns = ['CHR','SNP','BP']
snps_chr = snps_chr[['CHR','BP','SNP']]
snps_chr['CM'] = 0
snps_chr['base'] = 1
snps_chr = pd.concat([snps_chr, pd.DataFrame(index = snps_chr.index, columns = all_scores.columns, data = 0, dtype = float)], axis = 1)
ref_chr = ref.loc[ref.CHR == chrom,:]

idx = 0
tic = t()
for g, row in ref_chr.iterrows():
    if not g in all_scores.index: continue
    snps_chr.loc[(row.FROM <= snps_chr.BP) & (snps_chr.BP <= row.TO), all_scores.columns] += all_scores.loc[g,:]
    idx += 1
    print(f'{idx} / {ref_chr.shape[0]}, time = {t()-tic:.3f}')

snps_chr.to_csv(f'wang_2025_neocortex/{chrom}.annot', sep = '\t', index = False)
