import os
import sys

fid = sys.argv[1]

os.chdir('/rds/project/rds-Nl99R8pHODQ/ref')
files = [0]
f = open(f'freq_{fid}.vcf')

for chrom in range(1,25):
    files.append(open(f'freq_alfa_chr{chrom}.vcf','a'))

files.append(open('freq_alfa_insertions.vcf','a'))

line = '   '
while len(line) > 0:
    line = f.readline()
    try: tmp = line.split()
    except: continue
    
    line = line.replace('\n','')
    try:
        chrom = line[:9].replace('NC_','')
        chrom = int(chrom)
        print(line, file = files[chrom])
    except:
        print(line, file = files[-1])