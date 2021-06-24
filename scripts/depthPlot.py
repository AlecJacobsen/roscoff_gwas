import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", dest="input", help="Input File")
parser.add_argument("-s", "--samples", dest="samples", help="Samples File")
parser.add_argument("-o", "--output", dest="output", help="Output File")
args = parser.parse_args()

depth_file = args.input
samples_file = args.samples
out_file = args.output

depth = pd.read_csv(depth_file,sep = '\t')
sampleInfo = pd.read_csv(samples_file, header = None)
sampleInfo.columns = ['sample','barcode','INDV','plate','row','col']
depth = depth.merge(sampleInfo, on = ['INDV'])
depth = depth.sort_values(by = ['MEAN_DEPTH','INDV'])
col = []
plates = ['PL000093','PL000094','PL000095']
#plates = depth['plate'].unique() removed for consistancy
for line in depth['plate']:
    if line == plates[0]:
        col.append('firebrick')
    if line == plates[1]:
        col.append('teal')
    if line == plates[2]:
        col.append('darkorange')
plt.figure(figsize = (20,10))
plt.bar(depth['INDV'],depth['MEAN_DEPTH'], color = col)
legend_elements = [Patch(facecolor='firebrick',label=plates[0]),
                   Patch(facecolor='teal',label=plates[1]),
                   Patch(facecolor='darkorange',label=plates[2])]
plt.legend(handles=legend_elements,fontsize = 20, title = 'Plate', title_fontsize = 18)
plt.tick_params(axis='x', which='both',bottom=False,top=False,labelbottom=False)
plt.grid(axis='y', alpha = 0.25)
plt.title('Average Depth Per Individual', fontsize = 20)
plt.ylabel('Depth',fontsize = 18)
plt.xlabel('Individual',fontsize = 18)
plt.savefig(out_file)
