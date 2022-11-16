import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



import pandas as pd



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
bindingdb = bindingdb.drop(columns='Y')

bindingdb_negatives = bindingdb[bindingdb.Label == 0]
bindingdb_positives = bindingdb[bindingdb.Label == 1]

#bindingdb.groupby('UniprotID').sum()

grouped_bindingdb = bindingdb.groupby(by='UniprotID', group_keys=True).count()

selected_prot = grouped_bindingdb.sort_values(by='PubChemID', ascending=False).iloc[10:20].index.tolist()




protein = selected_prot[0]

bindingdb[bindingdb.UniprotID == protein]


possible_drugs_df = bindingdb[(bindingdb['UniprotID'] == protein) & (bindingdb['Label'] == 1)]

check_most_connect_drug = []
for drug in possible_drugs_df.PubChemID.tolist():
    total = bindingdb_positives[bindingdb_positives.PubChemID == drug].Label.sum()
    check_most_connect_drug.append((drug, total))

drug = max(check_most_connect_drug)[0]

print(f"Working with true pair {drug} - {protein}")

drop_pair = bindingdb[(bindingdb['UniprotID'] == protein) & (bindingdb['PubChemID'] == drug)].index.tolist()[0]


df_pairs = bindingdb[bindingdb.PubChemID == drug].drop(drop_pair)

# create a new column that gets the  rms value
test_ = df_pairs.UniprotID.tolist()[3] #in RMSD.index
test_
protein in RMSD.index

def get_rmsd_value(protein_selected, protein_check):
    if (protein_selected in RMSD.index) and (protein_check in RMSD.index):
        out = RMSD.loc[protein_selected].loc[protein_check]
    else:
        out = None
    return out

df_pairs['RMSD'] = df_pairs['UniprotID'].apply(lambda x: get_rmsd_value(protein, x))

df_pairs = df_pairs.dropna()

df_pairs

dict_col = {1: 'green', 0: 'black'}
colour = [dict_col.get(value) for value in df_pairs.Label.tolist()]
x_values = df_pairs.RMSD.tolist()
y_values = [75] * len(x_values)

#plt.scatter(x_values, y_values, c=colour, s =4)
# create new column for
#plt.savefig('test_hist_python.png')



protein


from collections import Counter


def vis_distr(drug, protein, RMSD, df_pairs):
    rmsd_distr = RMSD.loc[protein].drop(index=protein).values.tolist()
    #len(rmsd_distr)

    # plot stuff
    plt.clf()
    plt.title(f"{drug} - {protein} : 1")
    plt.hist(rmsd_distr, bins=140, color='skyblue')

    # vertical lines
    plt.axvline(x = 2.5, color = 'k')
    plt.axvline(x = 5, color = 'k')
    plt.axvline(x = 20, color = 'k')

    # props of points
    dict_col = {1: 'red', 0: 'black'}
    dict_alpha = {1: 1, 0: 0.25}
    dict_height = {1: 100, 0: 75}

    colour = [dict_col.get(value) for value in df_pairs.Label.tolist()]
    alphas = [dict_alpha.get(value) for value in df_pairs.Label.tolist()]

    x_values = df_pairs.RMSD.tolist()
    y_values = [dict_height.get(value) for value in df_pairs.Label.tolist()]#[75] * len(x_values)
    plt.text(65, 150, 'a')
    plt.scatter(x_values, y_values, c=colour, s =4, alpha=alphas)
    # create new column for
    plt.savefig('test_hist_python.png', dpi =330)

dict(Counter(df_pairs.Label.tolist()))

vis_distr(drug, protein, RMSD, df_pairs )

rmsd_vals = RMSD.values

va = rmsd_vals[np.triu_indices(rmsd_vals.shape[0] ,k = 1)]

np.triu(rmsd_vals, k=0)

np.triu_indices(rmsd_vals) 


def upper_tri_indexing(A):
    m = A.shape[0]
    r,c = np.triu_indices(m,1)
    return A[r,c]



vv = upper_tri_indexing(rmsd_vals)

plt.hist(vv)
plt.savefig('test_hist_python.png')