#!/usr/bin/env python3

'''
This script gets rsid from chromosome, base pair, and two alleles
Input:
  input text file, separated by whitespace and must have headers (!)
  column specification: 
    chromosome
    base pair (multiple possibilities!)
    ref allele, alt allele
    original ID
  genome build: hg19 or hg38
  reference file, use "%chr%" to put placeholder for each chromosome, "%build%" for genome build
  output prefix (auto appends '.txt' to the end)

Author:   Yuankai He (yh464@cam.ac.uk)
Date:     2025-08-10
'''
def num_chrom(df):
  if df.CHR.dtype != int: df['CHR'] = df.CHR.str.replace('chr','').replace('X','23').replace(
    'Y','24').replace('XY','25').replace('MT','26').replace('M','26').astype(int)
  return df

def load_ref_ucsc(file):
  import pandas as pd
  ref = pd.read_table(file, usecols = [1, 2, 3, 4, 7, 9], header = None)
  # col 11 = 'type'
  ref.columns = ['CHR','START','STOP','SNP','A1','A2']
  ref = num_chrom(ref)
  ref['A2'] = ref.A2.str.split('/')
  ref = ref.explode('A2')
  ref = ref.loc[(ref.A2 != 'lengthTooLong') & (ref.A2 != ref.A1),:]
  return ref, ['STOP','START']

def load_vcf(file):
  import pandas as pd
  ref = pd.read_table(file, usecols = [0,1,2,3,4], header = None, comment = '#')
  ref.columns = ['CHR','POS','SNP','A1','A2']
  ref = num_chrom(ref)
  ref['A2'] = ref.A2.str.split(',')
  ref = ref.explode('A2')
  ref = ref.loc[(ref.A2 != 'lengthTooLong') & (ref.A2 != ref.A1),:]
  return ref, ['POS']

def preprocess_indels(df):
  indel = df.A1.str.len() != df.A2.str.len()
  new_a1s = []; new_a2s = []
  for i in df.loc[indel, :].index:
    a1 = df.loc[i,'A1']; a2 = df.loc[i,'A2']
    if a1.startswith(a2): a1 = a1.removeprefix(a2); a2 = '-'
    elif a1.endswith(a2): a1 = a1.removesuffix(a2); a2 = '-'
    elif a2.startswith(a1): a2 = a2.removeprefix(a1); a1 = '-'
    elif a2.endswith(a1): a2 = a2.removesuffix(a1); a1 = '-'
    new_a1s.append(a1); new_a2s.append(a2)
  df.loc[indel, 'A1'] = new_a1s; df.loc[indel, 'A2'] = new_a2s

def merge_ids(df, ref, pos_cols = ['POS'], ref_cols = ['POS']):
  '''df should be normalised to have columns CHR, ID, A1, A2 and pos_cols
  CHR must be integers
  ref_cols and pos_cols should be sorted in order of preference'''
  import pandas as pd
  if not all([x in df.columns for x in pos_cols]): raise ValueError('pos_cols must be column names of df')
  if not isinstance(ref, pd.DataFrame): 
    if ref.find('.vcf') > 0: ref, ref_cols = load_vcf(ref); prep_indel = False
    else: ref, ref_cols = load_ref_ucsc(ref); prep_indel = True

  # retain a copy of the original data frame
  orig = df.copy()
  nsnp = df.shape[0]

  # UCSC format removes overlapping SNPs from alleles: e.g. A/AG -> -/G
  if prep_indel: df = preprocess_indels(df)

  # revert ref and alt alleles
  df_rev = df.copy()
  df_rev[['A1','A2']] = df[['A2','A1']]; df['INV'] = False; df_rev['INV'] = True
  df = pd.concat([df, df_rev], axis = 0); del df_rev

  # try different configs of genomic positions
  out_size = 0
  out = orig
  opt_pos_col = ''; opt_ref_col = ''
  for col in pos_cols:
    if out_size > 0.9 * nsnp: break
    for ref_col in ref_cols:
      tmp = pd.merge(df, ref, left_on = ['CHR',col,'A1','A2'], right_on = ['CHR',ref_col,'A1','A2'])
      if tmp.shape[0] > out_size: out = tmp; out_size = tmp.shape[0]; opt_pos_col = col; opt_ref_col = ref_col
      if out_size > 0.9 * nsnp: break
  del ref

  if out_size < 0.9 * nsnp:
    Warning(f'Only {out_size}/{nsnp} variants mapped, check genome build!\nReturning original IDs')
    out = orig; out['SNP'] = orig['ID']; return out, orig, None, None
  else:
    out = out[['CHR','ID','SNP']]
    out = pd.merge(orig, out, how = 'left') # report SNPs that are not matched
    missnp = out.loc[out.SNP.isna(),['CHR','ID','A1','A2'] + pos_cols] # snps that are not mapped
    out['SNP'] = out['SNP'].fillna(out['ID'])
    return out[['CHR','ID','A1','A2'] + pos_cols + ['SNP']], missnp, opt_pos_col, opt_ref_col
  
def main(args):
  import pandas as pd
  df = pd.read_table(args._in)

  # normalise column names
  df['ID'] = df[args.snp]
  df['CHR'] = df[args.chrom]
  df['A1'] = df[args.a1]; df['A2'] = df[args.a2]
  df = df[['ID','CHR','A1','A2'] + args.pos]
  df = num_chrom(df)

  # reference files
  ref = args.ref.replace(r'%build%',args.build)
  ref = [ref.replace(r'%chr%',str(x)) for x in range(1, 27)]

  # initialise
  opt_pos_col = None; opt_ref_col = None # update best-aligning coordinates for each chromosome
  out_match = []; out_mis = []; pos_cols = []; ref_cols = [] # and store in a list to test uniformity

  for chrom in df.CHR.unique():
    pos_col = args.pos; ref_col = ['STOP','START']
    if opt_pos_col != None: pos_col.remove(opt_pos_col); pos_col.insert(0, opt_pos_col)
    if opt_ref_col != None: ref_col.remove(opt_ref_col); ref_col.insert(0, opt_ref_col)

    df_chrom = df.loc[df.CHR == chrom,:]
    mat, mis, opt_pos_col, opt_ref_col = merge_ids(df_chrom, ref[chrom-1], pos_col, ref_col)
    out_match.append(mat); out_mis.append(mis); pos_cols.append(opt_pos_col); ref_cols.append(opt_ref_col)
  
  out_match = pd.concat(out_match); out_mis = pd.concat(out_mis)
  pos_cols = set(pos_cols); ref_cols = set(ref_cols)
  if len(pos_cols) > 1: Warning('Coordinates do not agree between chromosomes: ' + ' '.join(list(pos_cols)))
  if len(ref_cols) > 1: Warning('Different chromosomes align to different ref coordinates: ' + ' '.join(list(ref_cols)))
  print(f'Column {pos_cols[0]} of input file aligns to {ref_cols[0]} of the reference file')

  out_match.to_csv(f'{args.out}.txt', sep = '\t', index = False)
  out_mis.to_csv(f'{args.out}.mis.txt', sep = '\t', index = False)

if __name__ == '__main__':
  import argparse
  parser = argparse.ArgumentParser(description = 'This script fetches rsID from chr, pos and alleles')
  parser.add_argument('-i', '--in', dest = '_in', help = 'input file')
  parser.add_argument('-c', '--chrom', default = 'CHR', help = 'Chromosome column of input file')
  parser.add_argument('-p', '--pos', nargs = '+', help = 'Position columns of input files, multiple inputs')
  parser.add_argument('-a1', default = 'A1', help = 'Ref allele column of input file')
  parser.add_argument('-a2', default = 'A2', help = 'Alt allele column of input file')
  parser.add_argument('-s', '--snp', default = 'SNP', help = 'Original ID column of input file')
  parser.add_argument('-b','--build', choices = ['hg19','hg38'], default = 'hg38', help = 'Genome build')
  parser.add_argument('-r','--ref', help = 'Reference file, use %chr% as placeholder for chrom, %build% for build',
    default = '/rds/project/rds-Nl99R8pHODQ/ref/dbsnp/snp151_%build%_chr%chr%.txt.gz')
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