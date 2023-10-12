import matplotlib.pyplot as plt
import pandas as pd
import os

# Data
data = pd.read_csv(os.path.join('figs','datasets_statistics','dataset_stats.csv'), header=None)

# Bar width
bar_width = 0.35
index = range(len(datasets))

# Create a figure and axis
fig, ax = plt.subplots()

# Create bars for total nodes
bar1 = ax.bar(index, total_nodes, bar_width, label='Total Nodes', color='b')

# Create bars for total edges
bar2 = ax.bar([i + bar_width for i in index], total_edges, bar_width, label='Total Edges', color='g')

# Set the x-axis labels
ax.set_xticks([i + bar_width / 2 for i in index])
ax.set_xticklabels(datasets, rotation=45)

# Add labels and a legend
ax.set_xlabel('Dataset')
ax.set_ylabel('Count')
plt.title('Total Number of Nodes and Edges for Each Dataset')
plt.legend()

# Show the plot
plt.tight_layout()
plt.show()
