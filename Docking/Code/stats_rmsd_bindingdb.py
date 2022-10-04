
from turtle import pos
import pandas as pd

db_file_path = '../../DB/Data/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'


df = pd.read_csv(db_file_path, sep='\t')

bindingdb = pd.read_csv(db_file_path, sep="\t", header=0, usecols=['Drug_ID', 'Target_ID', 'Y'])
bindingdb = bindingdb.rename({'Drug_ID': 'PubChemID', 'Target_ID': 'UniprotID'}, axis=1)

bindingdb

bindingdb.loc[:, 'PubChemID'] = bindingdb.loc[:, 'PubChemID'].astype(int) # not float
bindingdb.loc[:, 'PubChemID'] = bindingdb.loc[:, 'PubChemID'].astype(str) # not float

bindingdb


threshold = 30
bindingdb['Label'] = [1 if x < threshold else 0 for x in bindingdb['Y']]
bindingdb = bindingdb.drop(columns='Y')

bindingdb_negatives = bindingdb[bindingdb.Label == 0]
bindingdb_positives = bindingdb[bindingdb.Label == 1]


RMSD_MATRIX = pd.read_pickle('../Data/pkls/RMSD_full_matrix.pkl')

def get_distr(dataset):
    grouped_biosnap = dataset.groupby(by='PubChemID', group_keys=True).count()
    no_multitarget = []
    two_targets = []
    multiple_targets = []
    i = 0
    for drug, n_proteins in grouped_biosnap.iterrows():
        i+=1
        print(f'{i} of {grouped_biosnap.shape[0]}', end='\r')
        # check if number of proteins != 1
        if n_proteins.values[0] == 1:
            aux = dataset[dataset.PubChemID == drug].UniprotID.tolist()[0]
            no_multitarget.append(aux)
        elif n_proteins.values[0] == 2:
            aux = dataset[dataset.PubChemID == drug].UniprotID.tolist()
            if aux[0] in RMSD_MATRIX.index and aux[1] in RMSD_MATRIX.index:
                rmsd_value = RMSD_MATRIX.loc[aux[0]][aux[1]]
                two_targets.append((aux, rmsd_value))
            else:
                two_targets.append((aux, None))
        else:
            lprots = dataset[dataset.PubChemID == drug].UniprotID.tolist()
            prot1 = lprots[0]
            rmsd_values = []
            for prot1 in lprots:
                vals = [RMSD_MATRIX.loc[prot1][prot] 
                        for prot in lprots 
                        if (prot != prot1) and 
                        (prot in RMSD_MATRIX.index) and(prot1 in RMSD_MATRIX.index)]
                rmsd_values.extend(vals)
                    #continue
            multiple_targets.append((drug,rmsd_values))
    return no_multitarget, two_targets, multiple_targets



no_multitarget_pos, two_targets_pos, multiple_targets_pos = get_distr(bindingdb_positives)
no_multitarget_neg, two_targets_neg, multiple_targets_neg = get_distr(bindingdb_negatives)



import matplotlib.pyplot as plt
from tqdm import tqdm

def visdistr(typ = 'pos'):
    if typ == 'pos':
        two_targets = two_targets_pos
        multiple_targets = multiple_targets_pos
    else:
        two_targets = two_targets_neg
        multiple_targets = multiple_targets_neg  
    plt.clf()
    plt.title(f'Dristribution RMSD per Pair Proteins {typ}')
    plt.hist([rmsd for _, rmsd in two_targets if rmsd != None], bins=100)
    plt.savefig(f'bindingdb_rmsds_pairs_{typ}.png')
    #
    plt.clf()
    plt.title(f'Dristribution RMSD per Multiple Proteins for drug: {typ}')
    for _, rmsd_values in tqdm(multiple_targets):
        _ = plt.hist(rmsd_values, bins=150, alpha=0.3)
    plt.savefig(f'bindingdb_rmsd_{typ}.png')


visdistr(typ='neg')



vals_neg = [i for _, values in multiple_targets_neg for i in values]
vals_pos = [i for _, values in multiple_targets_pos for i in values]

len(vals_neg), len(vals_pos)

plt.clf()
plt.hist(vals_neg, bins = 500, color = "skyblue", label='negatives')
plt.hist(vals_pos, bins = 500, color = "lightsalmon", label='positives')
#plt.yscale('log')
plt.legend()
plt.savefig('test_rmsd_bindings_2.png')


# import seaborn as sns
# plt.clf()
# fig, (ax1, ax2) = plt.subplots(1, 2)
# ax1 = sns.displot(vals_pos, kind="kde", fill=True, color =  "lightsalmon", rug=True)
# ax2 = sns.displot(vals_neg, kind="kde", fill=True, color =  "skyblue", rug=True)
# plt.savefig('test_rmsd_bindings_2.png')

