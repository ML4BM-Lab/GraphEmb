import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.patches import Patch



# data
df = pd.read_pickle('../Data/pkls/test_tani.pkl')
# annotations
df_annot = pd.read_pickle('../Data/pkls/df_annot.pkl')

df_annot=df_annot.fillna('NotCat')

##### histogram

plt.clf()
plt.title('hist drugbank (triang)')
keep = np.triu(np.ones(df.shape)).astype('bool').reshape(df.size)
vals = df.stack()[keep].to_numpy()
#sns.displot(vals, kind='kde')
sns.histplot(vals, kde=True)
plt.savefig('../Results/drug_test_hist.png',dpi=300)

# heatmap of full datafram
############################
## REEPS: 1
plt.clf()
plt.title('Tanimoto Similitude (DrugBank)')
sns.heatmap(df, xticklabels=False, yticklabels=False)
#sns.clustermap(df, metric='euclidian')
plt.savefig('../Results/drug_test_2.png',dpi=300)

# plt.imshow(np.array(df.values.tolist()).astype('float'))
# plt.imshow(df, cmap=plt.cm.get_cmap("Reds"), interpolation="nearest")
# plt.colorbar()

plt.clf()
plt.title('Tanimoto Similitude Cluster (DrugBank)')
#sns.heatmap(df, xticklabels=False, yticklabels=False)
pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)

sns.clustermap(df, metric='euclidean',  cmap = pal, xticklabels=False, yticklabels=False)

plt.savefig('../Results/drug_test_2_clust_nolab.png',dpi=400)


############################
## REEPS: 2 - Annotated

#########################
##### test clustermap

plt.clf()

superclass = df_annot.superclass
unique_superclass = superclass.unique().tolist()
colors_type = sns.color_palette("ch:s=.70,rot=-.70", 8)
lut = dict(zip(unique_superclass, colors_type))
row_colors = superclass.map(lut)
pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True) # general
sns.clustermap(df, cmap = pal,row_colors=row_colors, yticklabels=False, xticklabels=False)

plt.savefig('../Results/drug_test_annot_3.png',dpi=400)


####
plt.clf()
matrix = df
drug2superclass  = dict(zip(df_annot.PubChemID, df_annot.superclass))
matrix['superclass'] = matrix.index.map(drug2superclass)
matrix
#sns.heatmap(df, xticklabels=False, yticklabels=False)
#pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)
#sns.clustermap(df, metric='euclidean',  cmap = pal, xticklabels=False, yticklabels=False)
pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)
#dict_classes = {superclass_unique[i]:i+1 for i in range(len(df_annot.superclass.unique()))}
types = matrix['superclass']
sorted_types = types.unique().tolist()
sorted_types.sort()
#colors_type = sns.color_palette("Set2", len(sorted_types))
colors_type = sns.color_palette("tab20", len(sorted_types))
lut_type = dict(zip(sorted_types, colors_type))
row_colors = types.map(lut_type)
# labels = np.random.random_integers(0,5, size=50)
# lut = dict(zip(set(labels), sns.hls_palette(len(set(labels)), l=0.5, s=0.8)))
# row_colors = pd.DataFrame(labels)[0].map(lut)

#sns.clustermap(matrix.iloc[:,:-1],yticklabels=False,xticklabels=False)
sns.clustermap(matrix.iloc[:,:-1], cmap = pal,row_colors=row_colors, yticklabels=False, xticklabels=False)
handles = [Patch(facecolor=lut_type[name]) for name in lut_type]
plt.legend(handles, lut_type, title='superclasss', bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')
plt.savefig('../Results/drug_test_annot_3_drugbank_v2.png',dpi=400)


