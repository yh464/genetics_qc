#!/usr/bin/env python3
'''
Author: Yuankai He
Correspondence: yh464@cam.ac.uk
Version 1: 2025-01-14

Create allele reference file from 1000 genomes reference panel
'''

import pandas as pd
import os
os.chdir('/rds/project/rds-Nl99R8pHODQ/ref')
from argparse import ArgumentParser
import numpy as np

def parse_freq(entry):
    chunk = entry.split(':')
    try:
        total = int(chunk[0])
        alt = sum([int(y) for y in chunk[1].split(',')])
        return 1 - alt/total
    except:
        return np.nan

parser = ArgumentParser()
parser.add_argument('chrom', type = int)
args = parser.parse_args()

chrom = args.chrom
vcf = pd.read_table(f'1000g_ref_chr{chrom}.vcf.gz', comment = '#', header = None, usecols = [1,3,4])
ref = pd.read_table(f'snp151_hg19_chr{chrom}.txt.gz', header = None, usecols = [2,4,6]).drop_duplicates()
freq = pd.read_table(f'freq_alfa_chr{chrom}.vcf.gz', header = None, usecols = [1,3,4,9]) # hg19

ref.columns = ['POS','SNP','strand']
vcf.columns = ['POS','A1','A2']
freq.columns = ['POS','A1','A2','AF1']
freq['AF1'] = [parse_freq(x) for x in freq['AF1']]
ref['POS'] += 1 # important: dbsnp uses a different counting system
freq_rev = freq.copy()
freq_rev.columns = ['SNP','A2','A1','AF1']
freq_rev['AF1'] *= -1
freq_rev['AF1'] += 1
freq.drop('A2', inplace = True, axis = 1)
freq_rev.drop('A2', inplace = True, axis = 1)

out = pd.merge(vcf, ref)
out = pd.concat([pd.merge(out, freq), pd.merge(out, freq_rev)]).sort_values(by = 'POS')

out = out[['SNP','POS','A1','A2',
           'AF1',
           'strand']]
out.insert(loc = 0, column = 'CHR', value = chrom)
out.to_csv(f'dbsnp_1kg_alfa_ref_hg19_chr{chrom}.txt', sep = '\t', index = False)

hg38 = pd.read_table(f'snp151_hg38_chr{chrom}.txt.gz', usecols = [2,4]).drop_duplicates()
hg38.columns = ['POS','SNP']
out.drop('POS', axis = 1, inplace = True)
out = pd.merge(out, hg38)
out = out[['CHR','SNP','POS','A1','A2',
           'AF1',
           'strand']].sort_values(by = 'POS')
out.to_csv(f'dbsnp_1kg_alfa_ref_hg38_chr{chrom}.txt', sep = '\t', index = False)