import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt

#load dataframe
df = pd.read_csv(os.path.join('figs','time_consumption','time_consumption_benchmarking.csv'))
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
# annotdf = df.apply(lambda x: x.apply(lambda y: str(round(np.log2(y),2)) if isinstance(y, int) else y))
# valuesdf = df.apply(lambda x: x.apply(lambda y: round(np.log2(y),2) if isinstance(y, int) else 0))
    
df = df.replace('OOT',float('nan'))
df = df.replace('OOR',float('nan'))
df = np.log2(df)

# Create a figure
colors = ["#033aa5", "#bb6f7c", "#4970e3", "#b4b9e3", "#e6b5a9", "#d34c6a", "#8dd5a3", "#f0c78d"]
colors = ["#ca586f", "#70a845", "#8761cc", "#b49041", "#688bcc", "#cd5d39", "#4aac8d", "#c361aa"]

plt.figure(figsize=(10, 6))

# Plotting points and connecting lines
datasets = df.columns
for i,model in enumerate(df.index.tolist()):
    plt.plot(datasets, df.iloc[i,:].values, marker='o', label=model, color=colors[i])

# Adding labels and legend
plt.xlabel('Datasets')
plt.ylabel('Seconds (log2)')
plt.title('Time taken by different models on various datasets')
plt.legend()
plt.xticks(rotation=45, ha='right')  # Rotate x-axis labels for better readability

# Show the plot
plt.tight_layout()
plt.show()
plt.savefig(os.path.join('figs','time_consumption','time_consumption_figure.pdf'))