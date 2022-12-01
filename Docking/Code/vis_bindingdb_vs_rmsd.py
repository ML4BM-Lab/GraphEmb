import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from collections import Counter
import random


import pandas as pd

#bindingdb[bindingdb.UniprotID == 'P34976'].Label.sum() #shape

##### FIRST VERSION
# protein = selected_prot[1]
# bindingdb[bindingdb.UniprotID == protein] #shape
# #for protein in selected_prot:
# #bindingdb[bindingdb.UniprotID == protein]
# possible_drugs_df = bindingdb[(bindingdb['UniprotID'] == protein) & (bindingdb['Label'] == 1)]
# check_most_connect_drug = []
# for drug in possible_drugs_df.PubChemID.tolist():
#     total = bindingdb_positives[bindingdb_positives.PubChemID == drug].Label.sum()
#     check_most_connect_drug.append((drug, total))
# #drug = possible_drugs_df.PubChemID.tolist()[random.randint(0, possible_drugs_df.shape[0]-1)]
# drug = max(check_most_connect_drug)[0]
#####

def vis_distr(drug, protein, RMSD, df_pairs, mapcol = 'autumn'):
    rmsd_distr = RMSD.loc[protein].drop(index=protein).values.tolist()
    #len(rmsd_distr)

    # plot stuff
    plt.clf()
    #plt.ylim(0,450)
    plt.hist(rmsd_distr, bins=140, color='skyblue')

    # vertical lines
    plt.axvline(x = 2.5, color = 'k')
    plt.axvline(x = 5, color = 'k')
    plt.axvline(x = 20, color = 'k')

    # props of points
    dict_col = {1: 'red', 0: 'black'}
    dict_alpha = {1: 1, 0: 0.25}
    dict_height = {1: 100, 0: 75}

    num_items = dict(Counter(df_pairs.Label.tolist()))

    #colour = [dict_col.get(value) for value in df_pairs.Label.tolist()]
    colour = df_pairs.Y.tolist()
    colour = [col if col < 100 else 100 for col in colour]
    #colour =  [float(i)/sum(colour) for i in colour]
    alphas = [1] * len(colour) #[dict_alpha.get(value) for value in df_pairs.Label.tolist()]

    x_values = df_pairs.RMSD.tolist()
    [a,b,c,d] = get_intervals(x_values)
    y_values = [dict_height.get(value) for value in df_pairs.Label.tolist()]#[75] * len(x_values)
    #plt.text(65, 150, f'')
    plt.title(f"{drug} - {protein} : 1 \n 0: {num_items.get(0)} 1:{num_items.get(1)}\n {a:.2f}/{b:.2f}/{c:.2f}/{d:.2f}")
    plt.scatter(x_values, y_values, c=colour, s =4, alpha=alphas, cmap = mapcol)
    # create new column for
    plt.colorbar()
    plt.savefig(f'ttest.png', dpi =330)
    #plt.savefig(f'test_figures/hist_{protein}_{drug}.png', dpi =330)


# 
def get_intervals(x_values):
    a,b,c,d = 0,0,0,0
    tot = len(x_values)
    for val in x_values:
        if val >0 and val < 2.5:
            a+=1
        elif val > 2.5 and val < 5:
            b+=1
        elif val > 5 and val < 20:
            c+=1
        else:
            d+=1
    return np.array([a,b,c,d])/tot * 100

#x_values

def get_rmsd_value(protein_selected, protein_check):
    if (protein_selected in RMSD.index) and (protein_check in RMSD.index):
        out = RMSD.loc[protein_selected].loc[protein_check]
    else:
        out = None
    return out


############################
############################

RMSD = pd.read_pickle('../Results/RMSD_full_matrix.pkl')

db_file_path = '../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'

df = pd.read_csv(db_file_path, sep='\t')




bindingdb = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'Target_ID', 'Y'])
bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID'}, axis=1)

bindingdb

bindingdb.loc[:, 'PubChemID'] = bindingdb.loc[:, 'PubChemID'].astype(int) # not float
bindingdb.loc[:, 'PubChemID'] = bindingdb.loc[:, 'PubChemID'].astype(str) # not float

bindingdb

# filter proteins before the rest (proteins with available structure i nRMSD matrix)

threshold = 30
bindingdb['Label'] = [1 if x < threshold else 0 for x in bindingdb['Y']]
#bindingdb = bindingdb.drop(columns='Y')

bindingdb_negatives = bindingdb[bindingdb.Label == 0]
bindingdb_positives = bindingdb[bindingdb.Label == 1]

#bindingdb.groupby('UniprotID').sum()

#grouped_bindingdb_prot = bindingdb.groupby(by='UniprotID', group_keys=True).count()
grouped_bindingdb_drug = bindingdb_positives.groupby(by='PubChemID', group_keys=True).count()#.sum(numeric_only=True)


# select proteins with more drugs attached from 10th to 20th
#selected_prot = grouped_bindingdb_prot.sort_values(by='PubChemID', ascending=False).iloc[10:20].index.tolist()

selected_drugs = grouped_bindingdb_drug.sort_values(by='UniprotID', ascending=False).iloc[:30].index.tolist()
# try to select instead those with highest positives instead of whatever in general.
#selected_prot = grouped_bindingdb_label.sort_values(by='Label', ascending=False).iloc[10:20].index.tolist()

#drug = selected_drugs[0]
## SECOND VERSION WITH DRUGS


for drug in selected_drugs:
    
    # def generate_df_pairs(bindingdb, drug, RMSD):
    # if not protein: # and protein=None as argument
    pos_prot = bindingdb[(bindingdb.PubChemID == drug) & (bindingdb.Label == 1)] #.Label.sum()

    pos_prot_st = [prot for prot in pos_prot.UniprotID.tolist()  if prot in RMSD.index]

    protein = pos_prot_st[random.randint(0, len(pos_prot_st)-1)]
    
    if protein not in RMSD.index:
        continue 

    print(f"Working with true pair {drug} - {protein}")
    drop_pair = bindingdb[(bindingdb['UniprotID'] == protein) & (bindingdb['PubChemID'] == drug)].index.tolist()[0]
    df_pairs = bindingdb[bindingdb.PubChemID == drug].drop(drop_pair)

    df_pairs['RMSD'] = df_pairs['UniprotID'].apply(lambda x: get_rmsd_value(protein, x))
    df_pairs = df_pairs.dropna()
    #df_pairs.shape
    # return df_pairs
    vis_distr(drug, protein, RMSD, df_pairs, mapcol = 'cool')
    #sleep(2)




drug = '9829523'
protein = 'P07332'

vis_distr(drug, protein, RMSD, df_pairs, mapcol = 'gray')