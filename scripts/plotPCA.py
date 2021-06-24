import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

parser = argparse.ArgumentParser()
parser.add_argument("--vec", dest="vec", help="Eigenvec file")
parser.add_argument("--val", dest="val", help="Eigenval file")
parser.add_argument("--out", dest="out", help="Output File")
args = parser.parse_args()

vec_file = args.vec
val_file = args.val
out_file = args.out

vecs = pd.read_csv(vec_file, sep = ' ', header = None)
vals = pd.read_csv(val_file, header = None)[0].tolist()

ev1 = vecs[2].tolist()
ev2 = vecs[3].tolist()
colors = []
for indiv in vecs[1]:
    pop = indiv.split('_')[0]
    if pop == 'FM':
        colors.append('cornflowerblue')
    if pop == 'NM':
        colors.append('darkorange')

plt.figure(figsize = (10,10))
plt.scatter(ev1,ev2,c = colors, s = 25)
plt.title("PCA of Roscoff Individuals", fontsize = 20)
plt.xlabel('PC 1 ({:.2f}% variance)'.format(vals[0]), fontsize = 18)
plt.ylabel('PC 2 ({:.2f}% variance)'.format(vals[1]), fontsize = 18)
legend_elements = [Patch(facecolor='cornflowerblue',label='FM'),
                   Patch(facecolor='darkorange',label='NM')]
plt.legend(handles=legend_elements,fontsize = 16, title = 'Strain', title_fontsize = 16)
plt.savefig(out_file)
