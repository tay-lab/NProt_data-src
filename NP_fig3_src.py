# Import libraries
import numpy as np
import math
import random
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import ProxseqClasses as PC

# Directory path
myDir = "/path/to/data"
os.chdir(myDir)

#*****
mpl.rcdefaults()
# Set font to be arial
mpl.rc('font', **{'sans-serif':'Arial', 'size':12})
mpl.rcParams['mathtext.rm'] = 'sans' # to have non-italic greek letter, use r'$\mathrm{\alpha}$', does NOT work with f-string
mpl.rcParams['axes.titlesize'] = 12
#to store text as text, not as path
new_rc_params = {'text.usetex': False,
                 "svg.fonttype": 'none'}
mpl.rcParams.update(new_rc_params)
# Set default tick size
mpl.rcParams['xtick.major.size'] = 5.5
mpl.rcParams['ytick.major.size'] = 5.5
mpl.rcParams['xtick.minor.size'] = 2.5
mpl.rcParams['ytick.minor.size'] = 2.5
# Default legend settings
mpl.rcParams['legend.fancybox'] = False
mpl.rcParams['legend.edgecolor'] = 'k'

#%% Import real data
dge_real = pd.read_csv('example_dge.csv',index_col=0)

#plt.hist(dge_real.sum(axis=0))
# Filter cells
# Keep cells with > 10 UMIs and below 7k UMIs
dge_real = dge_real.loc[:,(dge_real.sum(axis=0)>10) & (dge_real.sum(axis=0)<7000)]
# Keep cells with more than 20 detected PLA products
dge_real = dge_real.loc[:,((dge_real>0).sum(axis=0)>20)]

#Create pla object 
pla_dge = PC.plaObject(dge_real)

#calculate protein complex from pla count
pla_dge.predictComplex()

#wide-format data preparation
T_1_express = pd.DataFrame({'frac':(pla_dge.complex_count>0).sum(axis=1)})/pla_dge.shape[1]
T_1_express['probeA'] = [s.split(':')[0] for s in T_1_express.index]
T_1_express['probeB'] = [s.split(':')[1] for s in T_1_express.index]
T_1_express = T_1_express.pivot(index='probeA', columns='probeB', values='frac')

my_avg_T = pd.DataFrame({'mean':(pla_dge.complex_count).mean(axis=1)})
my_avg_T['probeA'] = [s.split(':')[0] for s in my_avg_T.index]
my_avg_T['probeB'] = [s.split(':')[1] for s in my_avg_T.index]
my_avg_T = my_avg_T.pivot(index='probeA', columns='probeB', values='mean')

#Heatmap visulization Fig 3e&f
my_order = ['B7', 'CD147', 'CD28', 'CD3', 'CD4', 'CD45RA', 'HLADR', 'ICAM1',
            'LFA1', 'PD1', 'PDL1','IgG_28', 'IgG_29', 'IgG_30', ]
fig, ax = plt.subplots(figsize=(9,6), ncols=2)
sns.heatmap(T_1_express.loc[my_order,my_order], ax=ax[0],
            linewidths=1, square=True,
            cmap='viridis',
            vmax=1, cbar_kws={'label':"Fraction of expressing cells",
                              'ticks':np.arange(0,1.1,0.25),'shrink':0.4},)
sns.heatmap(my_avg_T.loc[my_order,my_order],square=True, linewidths=1, ax=ax[1],
            cmap=sns.cubehelix_palette(light=0.95, dark=0.15, as_cmap=True),
            vmin=0, vmax=100, cbar_kws={'label':"UMI per cell",'shrink':0.4})
for i in ax:
    i.set_xlabel("Probe B target")
    i.set_ylabel("Probe A target")
    i.tick_params(axis='both', length=0)
fig.tight_layout(h_pad=1.8)
fig.savefig(myDir+"fig 3e&f.svg", bbox_inches='tight', pad_inches=0)