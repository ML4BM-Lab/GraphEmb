import os
import pandas as pd

PATH = '/home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/'

all_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(os.path.join(PATH) ,DB)for f in filenames ]

for fl in files = 


a = pd.read_pickle('/home/margaret/data/jfuente/DTI/Input4Models/DTiGEMS/Data/DrugBank/Davis_Drug_SIDER_SideEffect.pickle')