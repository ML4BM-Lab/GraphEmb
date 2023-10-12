import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
import os
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


#load dataframe
df = pd.read_csv(os.path.join('figs','time_consumption_benchmarking.csv'))
df.index = df['Model/Dataset']
df.drop('Model/Dataset',axis=1, inplace=True)


# Function to convert time values to seconds
def time_to_seconds(time_str):
    if time_str == 'OOT' or time_str == 'OOR':
        return time_str
    parts = time_str.split()
    seconds = 0
    for part in parts:
        if 'd' in part:
            seconds += int(part.replace('d', '')) * 86400
        elif 'h' in part:
            seconds += int(part.replace('h', '')) * 3600
        elif 'm' in part:
            seconds += int(part.replace('m', '')) * 60
        else:
            seconds += int(part.replace('s', ''))
    return seconds

# Apply the time_to_seconds function to all columns
for col in df.columns:
    df[col] = df[col].apply(time_to_seconds)

# Create the heatmap
plt.figure(figsize=(7, 5))
annotdf = df.apply(lambda x: x.apply(lambda y: str(round(np.log2(y),2)) if isinstance(y, int) else y))
valuesdf = df.apply(lambda x: x.apply(lambda y: round(np.log2(y),2) if isinstance(y, int) else 0))
ax = sns.heatmap(valuesdf, annot=annotdf, cmap="Purples", fmt='', cbar_kws={'label': 'Log2 (seconds)'}, xticklabels=True, yticklabels=True, linewidths=1, linecolor='black')
# Change the size of x-label and y-label
plt.xticks(fontsize=4)  # Set the x-label font size to 12
plt.yticks(fontsize=6)  # Set the y-label font size to 14
#ax.xaxis.tick_top() 
plt.title("Time Consumption Heatmap (Log2 Scale)")
plt.show()
plt.savefig(os.path.join('figs','time_consumption_heatmap.pdf'))