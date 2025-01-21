#!/usr/bin/env python3
fin = open('psychencode_mixed_sqtl.txt')
fout = [0]
for chrom in range(1,26):
    fout.append(open(f'psychencode_sqtl/psychencode_sqtl_chr{chrom}_4besd.txt','a'))
fout.append(open('psychencode_sqtl/psychencode_sqtl_insertions_4besd.txt','a'))

line = '    '
columns = [0,1,2,-3,-2]
idx = 0
while len(line) > 0:
    line = fin.readline()
    idx += 1
    if idx % 10000000 == 0:
        print(f'{idx} lines read')
    
    try: 
        out = line.split()
        chrom = int(out[0].split(':')[0])
        out = '\t'.join([out[x] for x in columns])
        print(out, file = fout[chrom])
    except: 
        print(out, file = fout[-1])
