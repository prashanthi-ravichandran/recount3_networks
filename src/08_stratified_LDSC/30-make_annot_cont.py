import pandas as pd
import numpy as np
import argparse
import os
import glob
import sys
from intervaltree import Interval, IntervalTree

net = sys.argv[1]
centrality = sys.argv[2]
chr = sys.argv[3]
bedFile = "/data/abattle4/prashanthi/recount3_paper/results/s_LDSC/%s/%s/bed/%s.bed"%(net, centrality, net)
bimFile = "/data/abattle4/prashanthi/recount3_paper/data/s_LDSC/1000G_EUR_Phase3_plink/1000G.EUR.QC.%s.bim"%chr
annotfile_final = "/data/abattle4/prashanthi/recount3_paper/results/s_LDSC/%s/%s/ldscore/%s.%s.annot.gz"%(net, centrality, net, chr)
print(bedFile)
print(annotfile_final)
my_df = pd.read_csv(bedFile, sep='\t', names=['CHR','START','STOP','WEIGHT'])
my_df = my_df[my_df['CHR'] == 'chr%s' %chr]
my_df = my_df[['START','STOP','WEIGHT']]

df_bim = pd.read_csv(bimFile,delim_whitespace=True,usecols = [0,1,2,3],names = ['CHR','SNP','CM','BP'])

tree = IntervalTree.from_tuples(zip(
  my_df.START.values,
  my_df.STOP.values,
  my_df.drop(["START", "STOP"], axis=1).values.tolist() #save weight to list
))
annot = []
for row in df_bim.itertuples():
  ts = tree[row.BP]
if len(ts) == 0:
  annot.append(0)
else:
  max_val = 0
for it in ts:
  max_val = max(max_val, it.data[0])
annot.append(max_val)


df_bim['ANNOT'] = annot
#df_bim = df_bim.drop(["CHR","BP","SNP","CM"], axis = 1)
df_bim.to_csv(annotfile_final,sep='\t',compression='gzip',index=False)   
