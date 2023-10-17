import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from figs.datasets_statistics.scattermap import scattermap

def clean_plots():
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()

# Data
data = pd.read_csv(os.path.join('figs','datasets_statistics','dataset_stats.csv'), index_col=0).iloc[:-1,:]
demult = 10

# Create a figure and axis
fig, ax = plt.subplots(figsize=(8, 8))
ax = scattermap(data, linewidths=1, linecolor='black', marker_size=data/demult, cmap = 'flare')

# include legend
size_data = data.values
nsize = 4
intervals = (np.max(size_data) - np.min(size_data))/nsize
for i in range(nsize+1):
    realvalue = np.min(size_data) + i*intervals
    demultvalue = realvalue/demult
    ax.scatter(-1, -1, label=f"{realvalue:0.1f}", marker="o", c="r", s=demultvalue)
ax.legend(loc="upper left", prop={'size': 40})

# Show the plot
plt.tight_layout()
plt.grid(axis='y', linestyle='--', alpha=0.7)  # Add horizontal grid lines
plt.show()
plt.savefig(os.path.join('figs','datasets_statistics','dotplot.pdf'))
clean_plots()