#!/usr/bin/env python3
import sys
import os
import pandas as pd

prefix = sys.argv[1]

if not os.path.isdir(f'temp_snplist_{prefix}'): os.mkdir(f'temp_snplist_{prefix}')
if not os.path.isdir(prefix): os.mkdir(prefix)

snps = pd.read_table(f'{prefix}.esi', header = None)
for chrom in snps[0].unique():
    snps.loc[snps[0]==chrom, 1].to_csv(f'temp_snplist_{prefix}/chr{chrom}.txt', header = False, index = False)
    os.system(f'../../toolbox/smr --beqtl-summary {prefix} --extract-snp temp_snplist_{prefix}/chr{chrom}.txt --make-besd --out {prefix}/{prefix}_chr{chrom}')
    