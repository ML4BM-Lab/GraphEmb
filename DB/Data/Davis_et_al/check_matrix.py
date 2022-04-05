import pandas as pd
import numpy as np
import os

matrix = pd.read_csv(os.getcwd()+'/drug-target_interaction_affinities_Kd__Davis_et_al.2011.txt',sep=' ',header=None)
rows = pd.read_csv(os.getcwd()+'/drug_PubChem_CIDs.txt',header=None)
columns = pd.read_csv(os.getcwd()+'/target_gene_names.txt',header=None)
