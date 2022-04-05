import pandas as pd
import numpy as np
import os
from collections import Counter
import copy
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

#ff = pd.read_csv(os.getcwd()+'/BindingDB_All.tsv',sep='\t',error_bad_lines=False)
#From BindingDB_All, we filtered out those samples that did not have KD

##ONLY KD
ff_Kd =  pd.read_csv(os.getcwd()+'/BindingDB_Kd.tsv',sep='\t',index_col=0)

#First, lets filter out some columns
ff_Kd_filtered = ff_Kd[ff_Kd.columns[0:49]]

#FILTER BY HOMO SAPIENS
organisms = 'Target Source Organism According to Curator or DataSource'

p = []
for i,x in enumerate(ff_Kd_filtered.loc[:,organisms].values):
    if 'homo' in str(x).lower():
        p.append(i)
    #
    
#subset
ff_Kd_filtered_humans = ff_Kd_filtered.iloc[p,:]
ff_Kd_filtered_humans.to_csv(os.getcwd()+'/BindingDB_Kd_filtered_HomoSapiens.tsv',sep='\t')

#Visualize
ff_Kd_filtered_humans.head()

#Interesting columns
int_columns = ['BindingDB Target Chain  Sequence','Ligand SMILES','Kd (nM)','UniProt (SwissProt) Primary ID of Target Chain',
                'Curation/DataSource','DrugBank ID of Ligand','ChEMBL ID of Ligand','Link to Ligand in BindingDB','Link to Target in BindingDB', 'PubChem CID']
ff_Kd_filtered_columns_subset_humans = ff_Kd_filtered_humans[int_columns]


ff_Kd_filtered_columns_subset_humans.to_csv(os.getcwd()+'/BindingDB_Kd_filtered_columns_subset_HomoSapiens.tsv',sep='\t')

#EXAMPLE:
ff_Kd_filtered_columns_subset_humans = pd.read_csv(os.getcwd()+'/BindingDB_Kd_filtered_columns_subset_HomoSapiens.tsv',sep='\t')

#Counter
smiles_amino_c = Counter(list(zip(ff_Kd_filtered_columns_subset_humans['Ligand SMILES'],ff_Kd_filtered_columns_subset_humans['BindingDB Target Chain  Sequence'])))
print(sorted(smiles_amino_c.values()))

#We have seen that there are some pairs that have multiple Kd, so we are going to get the median of the Kds for those, just to keep 1 instead of having duplicities.
smiles_amino_d = {}
for i in ff_Kd_filtered_columns_subset_humans.index.values:
    #
    #generate the pair
    pair_l = '||'.join(list(ff_Kd_filtered_columns_subset_humans.loc[i,['Ligand SMILES','BindingDB Target Chain  Sequence']].values))
    #
    if not pair_l in smiles_amino_d:
        smiles_amino_d[pair_l] = [i]
    else:
        smiles_amino_d[pair_l] += [i]


def convert_str_to_float_and_sum(iter):
    #
    float_iter = []
    for value in iter:
        if value[0] == '>' or value[0] == '<':
            float_iter.append(float(value[1:]))
        else:
            float_iter.append(float(value))
    #
    return sum(float_iter)


#iterate
ff_Kd_filtered_columns_subset_humans_noduplicities = pd.DataFrame(columns=ff_Kd_filtered_columns_subset_humans.columns)
for i,pair in enumerate(smiles_amino_d):
    if len(smiles_amino_d[pair])>15:
        print(i)
        print(pair)
    #
    #append row
    new_row = copy.deepcopy(ff_Kd_filtered_columns_subset_humans.loc[smiles_amino_d[pair][0],int_columns])
    ff_Kd_filtered_columns_subset_humans_noduplicities.loc[i,:] = new_row
    #avg Kd and modify
    ff_Kd_filtered_columns_subset_humans_noduplicities.loc[i,'Kd (nM)'] = convert_str_to_float_and_sum(ff_Kd_filtered_columns_subset_humans.loc[smiles_amino_d[pair],'Kd (nM)'].values)

ff_Kd_filtered_columns_subset_humans_noduplicities.to_csv(os.getcwd()+'/BindingDB_Kd_filtered_columns_subset_HomoSapiens_noduplicities.tsv',sep='\t')

#print unique targets and drugs

print(f"Unique SMILES {len(set(ff_Kd_filtered_columns_subset_humans_noduplicities['Ligand SMILES'].values))}")
print(f"Unique Targets {len(set(ff_Kd_filtered_columns_subset_humans_noduplicities['BindingDB Target Chain  Sequence'].values))}")