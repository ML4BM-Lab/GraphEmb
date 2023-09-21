#/home/sevastopol/anaconda3/bin/python
import sys
import os
import subprocess
from Protein_kernels_SW import DAVIS_SW,  BINDINGDB_SW, DRUGBANK_SW, YAMANASHI_SW, BIOSNAP_SW
from Drugs_kernels_SIDER import DAVIS_SIDER,  BINDINGDB_SIDER, DRUGBANK_SIDER, YAMANASHI_SIDER, BIOSNAP_SIDER

## LETS COMPUTE THE "Pairwise matrix built with the Morgan Fingerprints with radius 2 (RDKIT)"
try:
    arg1 = sys.argv[1]
    try: 
        arg2 = sys.argv[2]
    except:
        arg2 = None
except:
    arg1, arg2 = None,None 

def get_paths(arg1,arg2):

    if arg2 == 'Yamanashi_et_al_GoldStandard':
        savepath = os.path.join(os.getcwd(),f'Methods/NeoDTI/Data/{arg2}/{arg1}/{arg1}_drug_')
    else:
        savepath = os.path.join(os.getcwd(),f'Methods/NeoDTI/Data/{arg1}/{arg1}_drug_')

    molpath = os.path.join(f'/tmp/{arg1}_drug_mol')
    protpath = os.path.join(f'/tmp/{arg1}_gene_symbol')
    dataset = arg2
    subdataset = arg1

    return molpath, protpath, savepath, dataset, subdataset

#get paths
molpath, protpath, savepath, dataset, subdataset = get_paths(arg1,arg2)

#-------------------------------------------------------------------------------------------------------------------------#
##DRUGS
#Generate MOL files
def GenMOLfiles(dataset,subdataset,path_to_mol):
    print('Generating MOL files for DRUGS')
    path_to_script = os.getcwd()+'/Methods/NeoDTI/Code/Preprocessing_genMOLfiles.py'
    #generate molfiles
    subprocess.run(["/home/sevastopol/anaconda3/bin/python",path_to_script, dataset, subdataset, path_to_mol],stdout=sys.stdout)

#MORGAN FINGERPRINT
def MorganFingerprint(molpath,savepath):
    print('Generating MorganFingerprint for DRUGS')
    path_to_script = os.getcwd()+'/Methods/NeoDTI/Code/Drugs_kernels_MorganFingerprint.py'
    subprocess.run(["/home/sevastopol/anaconda3/bin/python",path_to_script, molpath, savepath],stdout=sys.stdout)

#SIDER 
def SIDER(name):

    print(f'Generating SIDER for {name}')

    if name == 'Davis_et_al':
        DAVIS_SIDER()
    elif name == 'BIOSNAP':
        BIOSNAP_SIDER()
    elif name == 'DrugBank':
        DRUGBANK_SIDER()
    elif name == 'BindingDB':
        BINDINGDB_SIDER()
    elif name == 'Yamanashi_et_al_GoldStandard':
        YAMANASHI_SIDER()


#-------------------------------------------------------------------------------------------------------------------------#
##PROTEINS
#SW
def SW(name):

    print(f'Generating Smith Waterman for {name}')

    if name == 'Davis_et_al':
        DAVIS_SW()
    elif name == 'BIOSNAP':
        BIOSNAP_SW()
    elif name == 'DrugBank':
        DRUGBANK_SW()
    elif name == 'BindingDB':
        BINDINGDB_SW()
    elif name == 'Yamanashi_et_al_GoldStandard':
        YAMANASHI_SW()


#----------------------------------------------------------------------------------------------------------------------------#
#MAIN
#GenMOLfiles(dataset,subdataset,path_to_mol)
#MorganFingerprint(molpath,savepath)
