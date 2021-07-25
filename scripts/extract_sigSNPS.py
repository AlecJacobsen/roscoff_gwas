import argparse
import pandas as pd
import statsmodels.stats.multitest as multi

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", help="Input File")
args = parser.parse_args()

print('Reading in data')
data = pd.read_csv(args.input)

if 'P' and 'BP' in data.columns:
    data.rename(columns = {'P':'p', 'BP':'POS'}, inplace = True)
elif 'PVAL' in data.columns:
    data.rename(columns = {'PVAL':'p'}, inplace = True)

print('Extracting significant SNPs')
sig_snps = data.loc[multi.multipletests(data['p'])[0]]
sig_snps['CHR'] = "chr" + sig_snps['CHR'].astype(str)

print('Writing out results')
sig_snps[['CHR','POS']].to_csv(args.input.split('.')[0] + '.sig_snps', header = False, index = False, sep = '\t')

print('All done!')
