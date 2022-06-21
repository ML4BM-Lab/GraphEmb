import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


threshold = 5 # similarity

m = pd.read_csv('test_E.txt', header=0, sep ="\t+", index_col=0,  engine='python')
# m = pd.read_pickle('test_E_multiproces.pkl')
# read as pd.read_csv('test.txt', header=0, sep ="\t+", index_col=0)

plt.clf()
sns.heatmap(m<10, cmap='RdYlGn')
#m[m[m<5]>0.5]
plt.savefig('test_result_new.png')
