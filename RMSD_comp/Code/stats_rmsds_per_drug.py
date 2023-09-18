import numpy as np
import pandas as pd

import os
import pandas as pd
import numpy as np
import logging
import networkx as nx
import json
import helper_functions as hf


RMSD_MATRIX = pd.read_pickle('../Data/pkls/RMSD_full_matrix.pkl')



odtis = hf.original_dtis()

biosnap = odtis.biosnap()

grouped_biosnap = biosnap.groupby(by='DrugBankID', group_keys=True).count()

no_multitarget = []
two_targets = []
multiple_targets = []
i = 0
for drug, n_proteins in grouped_biosnap.iterrows():
    i+=1
    print(f'{i} of {grouped_biosnap.shape[0]}', end='\r')
    # check if number of proteins != 1
    if n_proteins.values[0] == 1:
        aux = biosnap[biosnap.DrugBankID == drug].UniprotID.tolist()[0]
        no_multitarget.append(aux)
    elif n_proteins.values[0] == 2:
        aux = biosnap[biosnap.DrugBankID == drug].UniprotID.tolist()
        if aux[0] in RMSD_MATRIX.index and aux[1] in RMSD_MATRIX.index:
            rmsd_value = RMSD_MATRIX.loc[aux[0]][aux[1]]
            two_targets.append((aux, rmsd_value))
        else:
            two_targets.append((aux, None))
    else:
        lprots = biosnap[biosnap.DrugBankID == drug].UniprotID.tolist()
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
    # create an else to also keep those 




import matplotlib.pyplot as plt
plt.clf()
plt.title('Dristribution RMSD per Pair Proteins')
plt.hist([rmsd for _, rmsd in two_targets if rmsd != None], bins=100)
plt.savefig('test_rmsds_pairs.png')


[(drug, rmsd) for drug, rmsd in two_targets if ( rmsd != None) and (rmsd > 20)]

np.mean([rmsd for _, rmsd in two_targets if rmsd != None])
np.std([rmsd for _, rmsd in two_targets if rmsd != None])


from tqdm import tqdm

plt.clf()
plt.title('Dristribution RMSD per Multiple Proteins for only one drug')
for _, rmsd_values in tqdm(multiple_targets):
    _ = plt.hist(rmsd_values, bins=150, alpha=0.3)
plt.savefig('test_rmsds.png')
