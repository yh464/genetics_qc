#!/usr/bin/env python3
'''
Subjects lists and genetic QC thresholds
'''

import os
import pandas as pd
df = pd.read_table(
    '/rds/project/rb643/rds-rb643-ukbiobank2/Data_Phenotype/DataFetch_20022024/ukb677594.tab',
    usecols = ['f.eid','f.22009.0.1','f.22009.0.2','f.21000.0.0']).set_index('f.eid')
df = df[['f.22009.0.1','f.22009.0.2','f.21000.0.0']]

# imaging samples
img = os.listdir('/rds/project/rb643/rds-rb643-ukbiobank2/Data_Imaging')
ls = []
for x in img:
    if x[:3] != 'UKB': continue
    if int(x[3:]) not in df.index: continue
    ls.append(int(x[3:]))

# ancestry
eur = [1,1001,1002,1003,1004]

summary = []

# overall mean/SD
def pc_std(df, name):
    mean_1 = df.iloc[:,0].mean(); std_1 = df.iloc[:,0].std()
    mean_2 = df.iloc[:,1].mean(); std_2 = df.iloc[:,1].std()
    nsub_5 = sum(
        (df.iloc[:,0] > mean_1 - 5*std_1) & (df.iloc[:,0] < mean_1 + 5*std_1) &\
        (df.iloc[:,1] > mean_2 - 5*std_2) & (df.iloc[:,1] < mean_2 + 5*std_2)
    )
    nsub_3 = sum(
        (df.iloc[:,0] > mean_1 - 3*std_1) & (df.iloc[:,0] < mean_1 + 3*std_1) &\
        (df.iloc[:,1] > mean_2 - 3*std_2) & (df.iloc[:,1] < mean_2 + 3*std_2)
    )
    ntotal = df.shape[0]
    return pd.DataFrame(dict(
        sub = [name], mean_1 = mean_1, std_1 = std_1, mean_2 = mean_2, std_2 = std_2,
        ntotal = ntotal, nsub_5 = nsub_5, nsub_3 = nsub_3
        ))

summary.append(pc_std(df, 'Whole UKB'))
df1 = df.loc[df.iloc[:,2].isin(eur),:]
summary.append(pc_std(df1, 'European ancestry'))
df1 = df.loc[ls,:]
summary.append(pc_std(df1, 'Imaging sample'))
df1 = df1.loc[df1.iloc[:,2].isin(eur), :]
summary.append(pc_std(df1, 'Imaging sample EUR'))
summary = pd.concat(summary)
print(summary)
summary.to_csv('/rds/project/rb643/rds-rb643-ukbiobank2/Data_Users/yh464/genqc/nsub.txt', sep = '\t', index = False)
