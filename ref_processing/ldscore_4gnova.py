# requires bitarray package

import os, tempfile
os.chdir(os.path.dirname(os.path.realpath(__file__)))
tmpdir = tempfile.mkdtemp()

for x in os.listdir():
    if not os.path.isdir(x): continue
    if x.startswith('.'): continue
    for y in range(1, 23):
        input_file = f'{x}/{y}.annot'
        
        # first format the annot files to remove the first 5 columns
        out_file = f'{tmpdir}/{x}_{y}.annot'

        cmd = "awk '{for (i=6; i<NF; i++) printf $i \"\\t\"; print $NF}'" + f" {input_file} > {out_file}"
        print(cmd)
        os.system(cmd)

    # then run GNOVA LDSC
    cmd = '/rds/project/rds-Nl99R8pHODQ/toolbox/gnova/bin/python /rds/project/rds-Nl99R8pHODQ/toolbox/gnova/GNOVA/gnova.py --bfile /rds/project/rds-Nl99R8pHODQ/ref/1000g_eur_ldsc/chr@ '+\
        f'--annot {tmpdir}/{x}_@.annot --save-ld {x}.gnova.ldscore --out /dev/null /rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/gcorr/ldsc_sumstats/disorders/adhd2025.sumstats '+\
        '/rds/project/rds-Q6dKROTNf6s/Data_Users/yh464/gcorr/ldsc_sumstats/structural_factors/cortical_expansion.sumstats'
    print(cmd)
    os.system(cmd)