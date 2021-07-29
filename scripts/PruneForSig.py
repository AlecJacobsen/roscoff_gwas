import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--results", dest="res", help="GWAS results")
parser.add_argument("-i", "--prune_in", dest="prune_in", help="PLink prune.in file")
parser.add_argument("-n", "--number", dest="n", help="Number of markers to keep", type = int)
parser.add_argument("-o", "--output", dest="out", help="Output file")
args = parser.parse_args()

n = args.n

with open(args.res,'r') as f:
    header = f.readline()

if header.__contains__('\t'):
    gwas_res = pd.read_csv(args.res, sep = '\t')
if header.__contains__(','):
    gwas_res = pd.read_csv(args.res, sep = ',')

if 'BP' in gwas_res.columns:
    gwas_res.rename(columns = {'BP':'POS'}, inplace = True)
if 'PVAL' in gwas_res.columns:
    gwas_res.rename(columns = {'PVAL':'P'}, inplace = True)
gwas_res.columns = [name.upper() for name in gwas_res.columns]

prune_in = pd.read_csv(args.prune_in, sep = ':', header = None)
prune_in.columns = ['CHR','POS']

print(len(prune_in['CHR']),"Markers read in")
print('Filtering for', n, 'most significant markers')
merged = gwas_res.merge(prune_in, how = 'inner', on=['CHR','POS'])
sigs = merged.nsmallest(n,'P')[['CHR','POS']].sort_values(['CHR','POS'], axis = 0)
print('Writing out results')
sigs.to_csv(args.out, sep = ':', header = False, index = False)
print('Done!')
