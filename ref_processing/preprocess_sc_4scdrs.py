import argparse
parser = argparse.ArgumentParser('Preprocessing single-cell data for scDRS')
parser.add_argument('-i', dest = '_in', help = 'Input h5ad file')
parser.add_argument('-o', dest = 'out', help = 'Output h5ad file')
parser.add_argument('-f', dest = 'force', help = 'Force overwrite', action = 'store_true')
args = parser.parse_args()

import os
dirname = os.path.dirname(args.out)
if not os.path.isdir(dirname): os.system(f'mkdir -p {dirname}')

import scanpy as sc
import scdrs

if not os.path.isfile(args.out) or args.force:
  adata = scdrs.util.load_h5ad(args._in)
  sc.pp.neighbors(adata, n_neighbors = 15, n_pcs = 20)
  scdrs.preprocess(adata)
  sc.write(args.out, adata)