## Generate all the adjacency matrices necessary for run DDR, for a specific dataset.
## Remind that we need:
    ## DRUGS
    ## For drugs, 8 different matrices, each from different methods:
        ## 5 from RchemCPP ✓
            # LAMBDA -> sd2gramSpectrum
            # MARGINALIZED -> sd2gram or sd2gramSpectrum
            # MINMAXTANIMOTO -> sd2gramSpectrum
            # SPECTRUM -> sd2gramSpectrum
            # TANIMOTO -> sd2gramSpectrum
        ## 2 from FDA DB (IN MARGARET)
            # AERSBIT -> Weighted cosine correlation coefficient
            # AERSFREQ -> Weighted cosine correlation coefficient
        ## 1 from SIDER (PENDING)
            # SIDER_METRIC -> Weighted cosine correlation coefficient
        ## 1 from SIMCOMP metric (IN MARGARET)

    ## PROTEINS
    ## For proteins, 6 different matrices, each from different methods:

        ## 2 from KeBABS R package ✓
            ## MISMATCH -> 4 levels
            ## SPECTRUM -> 2 levels
        ## 1 from BioMART DB package, from csbl.go R package 
        ## 1 from BioGRID DB, building a graph and computing the 'shortest hop distance' metric ✓
        ## 1 from Smith Waterman score (IN MARGARET)

import sys
import os
import subprocess
from Drugs_kernels_AERS_FDA import DAVIS_AERS_FDA, BINDINGDB_AERS_FDA, DRUGBANK_AERS_FDA, YAMANASHI_AERS_FDA, BIOSNAP_AERS_FDA
from Drugs_kernels_SIDER import DAVIS_SIDER,  BINDINGDB_SIDER, DRUGBANK_SIDER, YAMANASHI_SIDER, BIOSNAP_SIDER
from Drugs_kernels_SIMCOMP import DAVIS_SIMCOMP,  BINDINGDB_SIMCOMP, DRUGBANK_SIMCOMP, YAMANASHI_SIMCOMP, BIOSNAP_SIMCOMP
from Protein_kernels_SW import DAVIS_SW,  BINDINGDB_SW, DRUGBANK_SW, YAMANASHI_SW, BIOSNAP_SW

try:
    arg1 = sys.argv[1]
    try: 
        arg2 = sys.argv[2]
    except:
        arg2 = None
except:
    arg1, arg2 = None,None 


def get_paths(arg1,arg2):

    if arg1 == 'Yamanashi_et_al_GoldStandard':
        drug_savepath = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg2}/{arg2}_drug_')
        protein_savepath = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg2}/{arg2}_prot_')
        molpath = os.path.join(f'/tmp/{arg2}_drug_mol')
        protpath = os.path.join(f'/tmp/{arg2}_gene_symbol')
        fastapath = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg2}/{arg2}_Targets_AA_sequences.tsv')
    else:
        drug_savepath = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg1}_drug_')
        protein_savepath = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg1}_prot_')
        molpath = os.path.join(f'/tmp/{arg1}_drug_mol')
        protpath = os.path.join(f'/tmp/{arg1}_gene_symbol')
        fastapath = os.path.join(os.getcwd(),f'DDR/Data/{arg1}/{arg1}_Targets_AA_sequences.tsv')

    dataset = arg1
    subdataset = arg2

    return molpath, protpath, drug_savepath, protein_savepath, dataset, subdataset,fastapath

#get paths
molpath, protpath, drug_savepath, protein_savepath, dataset, subdataset, fastapath = get_paths(arg1,arg2)

##DRUGS
#----------------------------------------------------------------------------------------------------------#

#Preprocesing
#Generate MOL files
def GenMOLfiles(dataset,subdataset,path_to_mol):
    print('Generating MOL files for DRUGS')
    path_to_script = os.getcwd()+'/Methods/DDR/Code/Preprocessing_genMOLfiles.py'
    #generate molfiles
    subprocess.run(["/home/sevastopol/anaconda3/bin/python",path_to_script, dataset, subdataset, path_to_mol],stdout=sys.stdout)

#RchemCPP
def RchemCPP(path_to_mol,drug_savepath):
    path_to_script = os.getcwd()+'/Methods/DDR/Code/Drugs_kernels_Rchemcpp.R'
    #Rchemcpp
    subprocess.run(["Rscript",path_to_script, path_to_mol, drug_savepath],stdout=sys.stdout)

#FDA DB
def AERS_FDA(name,subdataset):

    print(f'Generating AERS FDA for {name}')
    if name == 'Davis_et_al':
        DAVIS_AERS_FDA()
    elif name == 'BIOSNAP':
        BIOSNAP_AERS_FDA()
    elif name == 'DrugBank':
        DRUGBANK_AERS_FDA()
    elif name == 'BindingDB':
        BINDINGDB_AERS_FDA()
    elif name == 'Yamanashi_et_al_GoldStandard':
        YAMANASHI_AERS_FDA(subdataset)

#SIDER 
def SIDER(name,subdataset):

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
        YAMANASHI_SIDER(subdataset)

#SIMCOMP
def SIMCOMP(name,subdataset):

    print(f'Generating SIMCOMP for {name}')

    if name == 'Davis_et_al':
        DAVIS_SIMCOMP()
    elif name == 'BIOSNAP':
        BIOSNAP_SIMCOMP()
    elif name == 'DrugBank':
        DRUGBANK_SIMCOMP()
    elif name == 'BindingDB':
        BINDINGDB_SIMCOMP()
    elif name == 'Yamanashi_et_al_GoldStandard':
        YAMANASHI_SIMCOMP(subdataset)

#----------------------------------------------------------------------------------------------------------#
##PROTEINS

#Preprocessing
#Generate Prot Symbols
def GenProtSymbols(dataset,subdataset):
    print('Generating Symbols for PROTEINS')
    path_to_script = os.getcwd()+'/Methods/DDR/Code/Preprocessing_genProtSymbols.py'
    #generate molfiles
    subprocess.run(["/home/sevastopol/anaconda3/bin/python",path_to_script, dataset, subdataset, protpath],stdout=sys.stdout)

#KeBABS
def KeBABS(fastapath,protein_savepath):
    #(NEED TARGET SEQUENCES)
    print(f'Generating KeBABS matrix')
    path_to_script = os.getcwd()+'/DDR/Code/Protein_kernels_kebabs.R'
    subprocess.run(["Rscript",path_to_script, fastapath,protein_savepath],stdout=sys.stdout)

#BioGRID
def BioGRID(protein_savepath,protpath):
    print(f'Generating BioGRID matrix')
    path_to_script = os.getcwd()+'/Methods/DDR/Code/Protein_kernels_BioGRID.py'
    subprocess.run(["/home/sevastopol/anaconda3/bin/python",path_to_script, protein_savepath, protpath],stdout=sys.stdout)

#BioMART
def BioMART(protpath,dataset,protein_savepath):
    #GenProtSymbols(name,dataset_name)
    print(f'Generating BioMART matrix')
    path_to_script = os.getcwd()+'/Methods/DDR/Code/Protein_kernels_bioMART.R'
    subprocess.run(["Rscript",path_to_script, protpath, dataset, protein_savepath],stdout=sys.stdout)

#SW
def SW(name,subdataset):

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
        YAMANASHI_SW(subdataset)

#---------------------------------------------------------------- MAIN ------------------------------------------#
##PREPARE
#GenMOLfiles(dataset,subdataset,molpath)
#GenProtSymbols(dataset,subdataset)

##Drugs
#RchemCPP(molpath,drug_savepath)
#AERS_FDA(dataset,subdataset)
#SIDER(dataset,subdataset)
#SIMCOMP(dataset,subdataset)

##Proteins
#BioMART(protpath,dataset,protein_savepath)
#BioGRID(protein_savepath,protpath)
#SW(dataset,subdataset)
#KeBABS(fastapath,protein_savepath)

