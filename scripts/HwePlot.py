import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", help="Input File")
parser.add_argument("-o", "--output", dest="output", help="Output File")
args = parser.parse_args()

hwe_file = args.input
out_file = args.output

hwe = pd.read_csv(hwe_file,sep = '\t')
plt.figure(figsize = (20,10))
plt.hist(hwe['P_HWE'],bins = 200, log = True)
plt.grid(axis='y', alpha = 0.25)
plt.title('Distribution of HWE for SNPs', fontsize = 20)
plt.ylabel('Proportion',fontsize = 18)
plt.xlabel('Signifigance (p)',fontsize = 18)
plt.savefig(out_file)
