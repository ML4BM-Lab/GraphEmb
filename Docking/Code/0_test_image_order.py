import pandas as pd
from tqdm import tqdm
import requests
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# threshold = 5 # similarity
# m = pd.read_csv('test_E.txt', header=0, sep ="\t+", index_col=0,  engine='python')
# # m = pd.read_pickle('test_E_multiproces.pkl')
# # read as pd.read_csv('test.txt', header=0, sep ="\t+", index_col=0)
# plt.clf()
# sns.heatmap(m<10, cmap='RdYlGn')
# #m[m[m<5]>0.5]
# plt.savefig('test_result_new.png')

matrix = pd.read_csv('test_E.txt', header=0, sep ="\t+", index_col=0,  engine='python')
# matrix = pd.read_pickle('test_E_multiproces.pkl')
matrix
#  matrix.applymap(float)
list_unis = matrix.columns.tolist()
list_unis

data_uniprot = pd.DataFrame(columns=['UniprotID', 'seq_len', 'keyword', 'EC'])
pop_words = ['3d-structure', 'calcium', 'zinc',
                'mental retardation', 'neurogenesis']
words = []
for uniprotid in tqdm(list_unis):
    #uniprotid
    r = requests.get(f'https://www.uniprot.org/uniprot/{uniprotid}.txt')
    a = r.text.split('\n')
    # SQ INFO
    SQ_info = None
    check_SQ = any('SQ' in string for string in a)
    SQ_info = [line.split() for line in a if 'SQ ' in line if check_SQ][0][2]
    ECode = None
    check_E = any('EC=' in string for string in a)
    check_BR = any('BRENDA;' in string for string in a)
    if check_E: 
        E_Code_line = [line.split() for line in a if 'EC=' in line][0]
        ECode = [element for element in E_Code_line if '=' in element]
        ECode = ECode[0].replace('EC=', '').strip(';')
        #Code
    elif check_BR:
        BR_Code_line = [line.split() for line in a if 'BRENDA;' in line]
        ECode = BR_Code_line[0][2].strip(';')
    # KEYWORDS
    check_KW = any('KW' in string for string in a)
    if check_KW:
        KW_info = [line.replace('KW', '').replace('.', '').split(';') for line in a if 'KW  ' in line]
        keywords = [x.strip().lower() for list_ in KW_info for x in list_]
        keywords = list(filter(lambda x: x != "", keywords)) # check non empty strings
        for key in keywords:
            if key in pop_words: keywords.remove(key)
            # special cases 
        if ECode:
            keyword = 'enzyme'
        elif 'ion channel' in keywords: 
            keyword = 'ion channel'
        elif 'g-protein coupled receptor' in keywords:
            keyword = 'gpcr'
        elif 'transport' in keywords:
            keyword = 'transport'
        elif 'chaperone' in keywords:
            keyword = 'chaperone'
        elif 'hormone' in keywords:
            keyword = 'hormone'
        elif 'transcription' in keywords:
            keyword = 'transcription'
        else:
            keyword = 'notAnnotated' #keywords
            words.append(keywords)
        #print('\n')
        #time.sleep(3)
    data_uniprot.loc[len(data_uniprot.index)] = [uniprotid, SQ_info, keyword, ECode]



data_uniprot

enz_dic = {
            1: 'Oxidoreductases',
            2: 'Transferases',
            3: 'Hydrolases',
            4: 'Lyases',
            5: 'Isomerases',
            6: 'Ligases',
            7: 'Translocases',
            8: 'notAnnotated'
}

data_uniprot[data_uniprot.keyword != 'enzyme']


data_uniprot.EC.fillna(value=np.nan, inplace=True)

#data_uniprot.EC = data_uniprot.EC.apply(lambda x: str(x)[0]).astype(float)
data_uniprot.EC  = data_uniprot.EC.dropna().apply(lambda x: str(x)[0]).astype(int)

data_uniprot = data_uniprot.sort_values(by='EC').reset_index()
ordered = data_uniprot.sort_values(by='EC').reset_index().UniprotID.tolist()

mat_or = pd.DataFrame(matrix, columns=ordered, index=ordered)

data_uniprot.EC.map(enz_dic)

import matplotlib.pyplot as plt
import seaborn as sns

mi_ = []
for type in range(1,8):
    mi = data_uniprot[data_uniprot.EC == type ].index.min()
    ma = data_uniprot[data_uniprot.EC == type ].index.max()
    (mi, ma)
    mi_.append(mi)

plt.clf()
# antes mat <5
# # mat_or = mat_or.applymap(float) 

ax = sns.heatmap(mat_or<5, annot=False,yticklabels=False,xticklabels=False, cmap='ocean')
#ax = sns.heatmap(matrix)

#mat_or[mat_or[mat_or<5]>0.5]
plt.text(5, 340, enz_dic[1], fontsize='x-small')
plt.text(95, 340, enz_dic[2], fontsize='x-small')
plt.text(190, 340, enz_dic[3], fontsize='x-small')

ax.hlines(mi_, *ax.get_xlim(), color='r', linestyles='dotted')
ax.vlines(mi_, *ax.get_ylim(), color='r', linestyles='dotted')

plt.savefig('test_result_ordered_2_new.png')

##### test clustermap
from matplotlib.patches import Patch

plt.clf()
pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)
types = matrix['type']
colors = sns.color_palette("ch:s=.70,rot=-.70", 8)

sorted_types = types.unique().tolist()
sorted_types.sort()
lut = dict(zip(sorted_types,colors ))
row_colors = types.map(lut)
#sns.clustermap(matrix.iloc[:,:-1],yticklabels=False,xticklabels=False)
sns.clustermap(matrix.iloc[:,:-1], cmap = pal, row_colors=row_colors, yticklabels=False, xticklabels=False)

handles = [Patch(facecolor=lut[name]) for name in lut]

plt.legend(handles, lut, title='Types',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

plt.savefig('test_cluster_t.png')


uni2ec = dict(zip(data_uniprot.UniprotID, data_uniprot.EC) )
matrix['type'] = matrix.index.map(uni2ec)
matrix.type = matrix.type.fillna(8)


uni2seqlen = dict(zip(data_uniprot.UniprotID, data_uniprot.seq_len) )
matrix['seqlen'] = matrix.index.map(uni2seqlen)

plt.clf()
plt.bar(matrix.index.tolist(), matrix['seqlen'])
plt.savefig('sequencehist.png')
# 

plt.clf()

import glob

list_pdbs = glob.glob('../../Data/Clean_from_PDB/*.pdb')
list_pdbs = [file.replace('../../Data/Clean_from_PDB/', '').replace('.pdb', '') for file in list_pdbs]
list_afold = glob.glob('../../Data/Clean_from_AFold/*.pdb')
list_afold = [file.replace('../../Data/Clean_from_AFold/', '').replace('.pdb', '') for file in list_afold]

list_pdbs

matrix['source'] = '-'
matrix.index.isin(list_afold)
matrix['source'][matrix.index.isin(list_afold)] = 'AFold'
matrix['source'][matrix.index.isin(list_pdbs)] = 'PDB'
matrix


plt.clf()
pal = sns.color_palette("ch:s=.25,rot=-.25", as_cmap=True)

types = matrix['type']
colors_type = sns.color_palette("ch:s=.70,rot=-.70", 8)
sorted_types = types.unique().tolist()
sorted_types.sort()
lut_type = dict(zip(sorted_types,colors_type ))
row_colors = types.map(lut_type)

sources = matrix['source']
colors_sources = ['#b796dc', '#5aa2f0']
sorted_sources = sources.unique().tolist()
lut_sorted = dict(zip(sorted_sources, colors_sources  ))
colum_colors = sources.map(lut_sorted)

#sns.clustermap(matrix.iloc[:,:-1],yticklabels=False,xticklabels=False)
sns.clustermap(matrix.iloc[:,:-2], cmap = pal, row_colors=row_colors, col_colors=colum_colors ,yticklabels=False, xticklabels=False)

handles = [Patch(facecolor=lut_type[name]) for name in lut_type]
plt.legend(handles, lut_type, title='Types',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

handles = [Patch(facecolor=lut_sorted[name]) for name in lut_sorted]
plt.legend(handles, lut_sorted, title='Source',
           bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure, loc='upper right')

plt.savefig('test_cluster_t2.png')
