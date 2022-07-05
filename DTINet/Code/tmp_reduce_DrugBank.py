import numpy as np
import pandas as pd

dtis = pd.read_csv("../Data/DrugBank_big/DrugBank_DTIs.tsv", sep="\t")
dtis.columns = ['DrugBankID', 'UniprotID']

drugs = dtis.DrugBankID.unique().tolist()
proteins = dtis.UniprotID.unique().tolist()

# would need to check if fasta & smile
# then we can read big txt files

big_drugs = np.loadtxt("../Data/DrugBank_big/drug.txt", dtype='str').tolist()
big_proteins = np.loadtxt("../Data/DrugBank_big/protein.txt", dtype='str').tolist()

list_drug_nodes = list(set(drugs).intersection(set(big_drugs)))
list_protein_nodes = list(set(proteins).intersection(set(big_proteins)))