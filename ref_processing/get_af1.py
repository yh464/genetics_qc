#!/usr/bin/env python3
'''
This script gets rsid from chromosome, base pair, and two alleles
Input:
  input text file, separated by whitespace and must have headers (!)
  column specification: 
    rsID/SNP
    ref allele, alt allele
  reference file
  output prefix (auto appends '.txt' to the end)

Author:   Yuankai He (yh464@cam.ac.uk)
Date:     2025-08-20
'''
def num_chrom(df):
  if df.CHR.dtype != int: df['CHR'] = df.CHR.str.replace('chr','').replace('X','23').replace(
    'Y','24').replace('XY','25').replace('MT','26').replace('M','26').astype(int)
  return df

def load_ref_alfa(file, ancestry = 'EUR'):
  match ancestry:
    case 'EUR': ancestry_col = 9
    case 'EAS': ancestry_col = 11
    case 'AMR': ancestry_col = 14
    case 'SAS': ancestry_col = 16
    case 'AFR': ancestry_col = 18
    case _: ancestry_col = 15
  
  import pandas as pd
  df = pd.read_table(file, comment = '#', header = None, usecols = [2, 3, ancestry_col])
  counts = df[ancestry_col].str.split(':', expand = True).iloc[:,:2]
  counts[0] = counts[0].astype(int)
  counts[1] = counts[1].str.split(',', expand = True).fillna('0').astype(int).sum(axis = 1)
  df.columns = ['SNP','REF','info']
  df['AF2'] = counts[1]/counts[0]
  df['AF1'] = 1 - df.AF2
  return df.drop('info', axis = 1)

def merge_freq(df, ref, ancestry = 'EUR'):
  import pandas as pd
  if not isinstance(ref, pd.DataFrame):
    ref = load_ref_alfa(ref, ancestry)
  out = pd.merge(df, ref, how = 'left', on = 'SNP')
  invert = (~out.REF.isna()) & (out.REF != ref.A1)
  out.loc[invert,'AF1'] = out.loc[invert,'AF2']
  return out.drop('AF2', 'REF')

def main(args):
  import pandas as pd
  import numpy as np
  df = pd.read_table(args._in)

  # normalise column names
  df = df.rename(columns = {args.snp:'SNP', args.chrom:'CHR', args.a1:'A1',args.a2:'A2',
    'AF1': 'AF1_orig', 'AF2': 'AF2_orig'})
  df = num_chrom(df)
  
  df_by_chr = [df.loc[df.CHR == x,:] for x in df.CHR.unique()]
  ref_by_chr = [args.ref.replace('%chr%',str(x)) for x in df.CHR.unique]

  out = pd.concat([merge_freq(a, b, args.ancestry) for a,b in zip(df_by_chr, ref_by_chr)])
  out.to_csv(f'{args.out}.txt', sep = '\t', index = False)
  out.loc[out.AF1.isna(),:].to_csv(f'{args.out}.mis.txt', sep = '\t', index = False)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description = 'This script fetches allele frequencies from SNP and alleles')
  parser.add_argument('-i', '--in', dest = '_in', help = 'input file')
  parser.add_argument('-c', '--chrom', default = 'CHR', help = 'Chromosome column of input file')
  parser.add_argument('-a1', default = 'A1', help = 'Ref allele column of input file')
  parser.add_argument('-a2', default = 'A2', help = 'Alt allele column of input file')
  parser.add_argument('-s', '--snp', default = 'SNP', help = 'Original ID column of input file')
  parser.add_argument('-e','--eth', help = 'Ethnic group', choices = ['EUR','EAS','AMR','SAS','AFR','MIX'])
  parser.add_argument('-r','--ref', help = 'Reference file, use %chr% as placeholder for chrom',
    default = '/rds/project/rds-Nl99R8pHODQ/ref/alfa/freq_alfa_chr%chr%.vcf.gz')
  parser.add_argument('-o', '--out', help = 'output prefix', required = True)
  args = parser.parse_args()

  if args.out.endswith('.txt'): args.out = args.out.removesuffix('.txt')
  import os
  for attr in ['_in','out']: setattr(args, attr, os.path.realpath(getattr(args, attr)))
  from _utils import cmdhistory, logger
  logger.splash(args)
  cmdhistory.log()
  try: main(args)
  except: cmdhistory.errlog()