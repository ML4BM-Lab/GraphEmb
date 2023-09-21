import sys
import os
import subprocess
import pandas as pd
import numpy as np
from tqdm import tqdm

def postprocessing_AERS(arg1,arg2):

    if arg1 == 'Yamanashi_et_al_GoldStandard':
        fpath1 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg2}/{arg2}_drug_FDA_aers_bit.tsv')
        fpath2 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg2}/{arg2}_drug_FDA_aers_freq.tsv')
        fsave1 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg2}/{arg2}_drug_adjmat_FDA_aers_bit.tsv')
        fsave2 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg2}/{arg2}_drug_adjmat_FDA_aers_freq.tsv')
    else:
        fpath1 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg1}_drug_FDA_aers_bit.tsv')
        fpath2 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg1}_drug_FDA_aers_freq.tsv')
        fsave1 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg1}_drug_adjmat_FDA_aers_bit.tsv')
        fsave2 = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg1}_drug_adjmat_FDA_aers_freq.tsv')

    def tanimoto_coef(a,b):

        seta = set(a)
        setb = set(b)

        try:
            return len(seta.intersection(setb)) / len(seta.union(setb))
        except ZeroDivisionError:
            return 0
    
    def apply_postpro(fpath,fsave):
        #read mat
        print("Reading matrix")
        f1 = pd.read_csv(fpath,sep='\t',index_col=0)
        #get drug names
        dnames = sorted(f1.index)
        dnames_L = len(dnames)

        #build dicts
        name_aers_d = {}
        for i in tqdm(range(f1.shape[0]),desc='Building dict'):
            dname = f1.iloc[[i],:].index[0]
            aerv = [i for i,x in enumerate(f1.iloc[i,:].values) if x]
            #add to dict
            name_aers_d[dname] = aerv

        #init mat
        drug_drug_aers = np.ones(shape=(dnames_L,dnames_L))

        for i in tqdm(range(dnames_L),desc='Computing adjmat'):
            for j in range(i+1,dnames_L):
                drug_drug_aers[i,j] = tanimoto_coef(name_aers_d[dnames[i]],name_aers_d[dnames[j]])

        #get indices
        utriindex = np.triu_indices(drug_drug_aers.shape[0],k=1)

        #lowertri == uppertri
        drug_drug_aers.T[utriindex] = drug_drug_aers[utriindex]

        #build df
        aers_adjmat_df = pd.DataFrame(drug_drug_aers,index=dnames,columns=dnames)
    
        #save 
        aers_adjmat_df.to_csv(fsave,sep='\t')

    #aers bit
    print(f'Computing aers bit adjmat for {arg1} {arg2}')
    apply_postpro(fpath1,fsave1)
    #aers freq
    print(f'Computing aers freq adjmat for {arg1} {arg2}')
    apply_postpro(fpath2,fsave2)

try:
    arg1 = sys.argv[1]
    try: 
        arg2 = sys.argv[2]
    except:
        arg2 = ''
except:
    arg1, arg2 = '','' 

postprocessing_AERS(arg1,arg2)
