import os
import pandas as pd
import requests as r
import xml.etree.ElementTree as ET
from tqdm import tqdm
import numpy as np
import sys
from rdkit import Chem
from tqdm import tqdm

def getMol_DrugBank(drugID):
    req = r.get('https://go.drugbank.com/structures/small_molecule_drugs/'+drugID+'.mol')

    if req.status_code != 200 or 'END' not in req.text:
        print('.mol file for Drug '+ drugID+' not found')
        return ''
    else:
        return req.text+'\n'

def getMol_Kegg(drugID):
    #
    req = r.get('https://www.genome.jp/entry/-f+m+'+drugID)
    #
    if req.status_code != 200 or 'END' not in req.text:
        print('.mol file for Drug '+ drugID+' not found')
        return ''
    else:
        return req.text+'\n'

def generate_MOL_from_SMILE(pubchemid,smile,path_to_mol):
    path_to_file = path_to_mol+'/'+str(pubchemid)+'.mol'
    if not os.path.exists(path_to_file):
        with open(path_to_mol+'/'+str(pubchemid)+'.mol','w') as w:
            w.write(Chem.MolToMolBlock(Chem.MolFromSmiles(smile)))

##YAMANISHI
def genMOLfiles(dfname='Yamanashi_et_al_GoldStandard',name='NR',forcesave=False):
    ###---------------------------------------- Define paths ------------------------------------###
    #DTI
    if dfname == 'Yamanashi_et_al_GoldStandard':
        path = os.getcwd()+'/DBs/'+dfname+'/'+name+'/'+name.lower()+'_admat_dgc_mat_2_line.txt'
        path_to_mol = '/tmp/'+name+'_drug_mol'
    elif dfname == 'BindingDB':
        path = os.getcwd()+'/DBs/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.csv'
        path_to_mol = '/tmp/BindingDB_drug_mol'
    elif dfname == 'Davis_et_al':
        path = os.getcwd()+'/DBs/'+dfname+'/tdc_package_preprocessing/'+'DAVIS_et_al.tsv'
        path_to_mol = '/tmp/'+dfname+'_drug_mol'
    elif dfname == 'BIOSNAP':
        path = os.getcwd()+'/DBs/'+dfname+'/ChG-Miner_miner-chem-gene/'+'ChG-Miner_miner-chem-gene.tsv'
        path_to_mol = '/tmp/'+dfname+'_drug_mol'
    elif dfname == 'DrugBank':
        path = os.getcwd()+'/DBs/'+dfname+'/All_DB_Did_Tid_HumanSingleProteins_FDA_D_withSMI_T_NoFalseAA_format.tsv'
        path_to_mol = '/tmp/'+dfname+'_drug_mol'
    #
    #Load the dataset and obtain the Drug

    #Create empty folder
    if not os.path.isdir(path_to_mol):
        os.mkdir(path_to_mol)
    
    if dfname == 'Yamanashi_et_al_GoldStandard':

        drugs =  np.unique(pd.read_csv(path,sep='\t',header=None,names=['Drug','Protein'])['Drug'].values)
        #Generate a .mol file for each folder
        for drug in tqdm(drugs):
            #
            save_path = path_to_mol+'/'+drug+'.mol'
            if not os.path.isfile(save_path) or forcesave:
                drug_mol = getMol_Kegg(drug)
                #
                if len(drug_mol):
                    with open(save_path, "w") as f:
                        f.writelines(drug_mol)

    elif dfname == 'BindingDB':
        df = pd.read_csv(path,sep='\t',index_col=0)
        drugs_pubchem = sorted(set(map(int,df['Drug_ID'].values)))
        pubchem_smiles_dict = dict(zip(list(map(int,df.Drug_ID)),df.SMILES))
        #generate MOL from SMILE
        for pubchem_id in tqdm(drugs_pubchem):
            generate_MOL_from_SMILE(pubchem_id,pubchem_smiles_dict[pubchem_id],path_to_mol)
            
    elif dfname == 'Davis_et_al': #UNFINISHED
        df = pd.read_csv(path,sep='\t',index_col=0)
        drugs_pubchem = sorted(set(map(int,df.Drug_ID.values)))
        pubchem_smiles_dict = dict(zip(list(map(int,df.Drug_ID)),df.SMILES))
        #generate MOL from SMILE
        for pubchem_id in tqdm(drugs_pubchem):
            generate_MOL_from_SMILE(pubchem_id,pubchem_smiles_dict[pubchem_id],path_to_mol)

    elif dfname == 'BIOSNAP' or dfname == 'DrugBank':
        drugs =  sorted(np.unique(pd.read_csv(path,delimiter='\t',index_col=0).Drug.values))
        #
        for drug in tqdm(drugs):
            save_path = path_to_mol+'/'+drug+'.mol'
            if not os.path.isfile(save_path) or forcesave:
                drug_mol = getMol_DrugBank(drug)
                #
                if len(drug_mol):
                    with open(save_path, "w") as f:
                        f.writelines(drug_mol)

#generate mol files for name and dfname
genMOLfiles(sys.argv[1],sys.argv[2])
import os
import pandas as pd
import requests as r
import xml.etree.ElementTree as ET
from tqdm import tqdm
import numpy as np
import sys
from rdkit import Chem
from tqdm import tqdm

def getMol_DrugBank(drugID):
    req = r.get('https://go.drugbank.com/structures/small_molecule_drugs/'+drugID+'.mol')

    if req.status_code != 200 or 'END' not in req.text:
        print('.mol file for Drug '+ drugID+' not found')
        return ''
    else:
        return req.text+'\n'

def getMol_Kegg(drugID):
    #
    req = r.get('https://www.genome.jp/entry/-f+m+'+drugID)
    #
    if req.status_code != 200 or 'END' not in req.text:
        print('.mol file for Drug '+ drugID+' not found')
        return ''
    else:
        return req.text+'\n'

def generate_MOL_from_SMILE(pubchemid,smile,path_to_mol):
    path_to_file = path_to_mol+'/'+str(pubchemid)+'.mol'
    if not os.path.exists(path_to_file):
        with open(path_to_mol+'/'+str(pubchemid)+'.mol','w') as w:
            w.write(Chem.MolToMolBlock(Chem.MolFromSmiles(smile)))

##YAMANISHI
def genMOLfiles(dataset = 'Yamanashi_et_al_GoldStandard', subdataset = 'NR', path_to_mol = '/tmp/', forcesave = False):
    ###---------------------------------------- Define paths ------------------------------------###
    #DTI
    if dataset == 'Yamanashi_et_al_GoldStandard':
        path = os.getcwd()+'/DBs/'+dataset+'/'+subdataset+'/'+subdataset.lower()+'_admat_dgc_mat_2_line.txt'
    elif dataset == 'BindingDB':
        path = os.getcwd()+'/DBs/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
    elif dataset == 'Davis_et_al':
        path = os.getcwd()+'/DBs/'+dataset+'/tdc_package_preprocessing/'+'DAVIS_et_al.tsv'
    elif dataset == 'BIOSNAP':
        path = os.getcwd()+'/DBs/'+dataset+'/ChG-Miner_miner-chem-gene/'+'ChG-Miner_miner-chem-gene.tsv'
    elif dataset == 'DrugBank':
        path = os.getcwd()+'/DBs/'+dataset+'/DrugBank_DTIs.tsv'
    
    #Load the dataset and obtain the Drug

    #Create empty folder
    if not os.path.isdir(path_to_mol):
        os.mkdir(path_to_mol)
    
    if dataset == 'Yamanashi_et_al_GoldStandard':
        drugs =  np.unique(pd.read_csv(path,sep='\t',header=None,names=['Drug','Protein'])['Drug'].values)
        
        #Generate a .mol file for each folder
        for drug in tqdm(drugs):
            #
            save_path = path_to_mol+'/'+drug+'.mol'
            if not os.path.isfile(save_path) or forcesave:
                drug_mol = getMol_Kegg(drug)
                #
                if len(drug_mol):
                    with open(save_path, "w") as f:
                        f.writelines(drug_mol)

    elif dataset == 'BindingDB':
        df = pd.read_csv(path,sep='\t',index_col=0)
        drugs_pubchem = sorted(set(map(int,df['Drug_ID'].values)))
        pubchem_smiles_dict = dict(zip(list(map(int,df.Drug_ID)),df.SMILES))
        #generate MOL from SMILE
        for pubchem_id in tqdm(drugs_pubchem):
            generate_MOL_from_SMILE(pubchem_id,pubchem_smiles_dict[pubchem_id],path_to_mol)
            
    elif dataset == 'Davis_et_al': #UNFINISHED
        df = pd.read_csv(path,sep='\t',index_col=0)
        drugs_pubchem = sorted(set(map(int,df.Drug_ID.values)))
        pubchem_smiles_dict = dict(zip(list(map(int,df.Drug_ID)),df.SMILES))
        #generate MOL from SMILE
        for pubchem_id in tqdm(drugs_pubchem):
            generate_MOL_from_SMILE(pubchem_id,pubchem_smiles_dict[pubchem_id],path_to_mol)

    elif dataset == 'BIOSNAP' or dataset == 'DrugBank':
        drugs =  sorted(np.unique(pd.read_csv(path,delimiter='\t',index_col=0).Drug.values))
        #
        for drug in tqdm(drugs):
            save_path = path_to_mol+'/'+drug+'.mol'
            if not os.path.isfile(save_path) or forcesave:
                drug_mol = getMol_DrugBank(drug)
                #
                if len(drug_mol):
                    with open(save_path, "w") as f:
                        f.writelines(drug_mol)

#generate mol files for subdataset and dataset
dataset, subdataset, path_to_mol = sys.argv[1], sys.argv[2], sys.argv[3]
genMOLfiles(dataset, subdataset, path_to_mol)
