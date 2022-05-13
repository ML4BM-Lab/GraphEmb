import numpy as np
from tqdm import tqdm
import pandas as pd
import requests as r
import os
import sys

def getSymbol_Kegg(geneID):
        #
        req = r.get('http://rest.kegg.jp/list/'+geneID)
        #
        if req.status_code != 200 or 'hsa' not in req.text:
            print('symbol file for gene  '+ geneID+' not found')
            return ''
        else:
            return req.text.split('\t')[1].split(';')[0].split(', ')

def getSymbol(subdataset='E',dataset='Yamanashi_et_al_GoldStandard',path_to_gsymbol = '/tmp/'):

    #get path
    if dataset == 'Yamanashi_et_al_GoldStandard':
        path = os.getcwd()+'/DBs/'+dataset+'/'+subdataset+'/'+subdataset.lower()+'_admat_dgc_mat_2_line.txt'
    elif dataset == 'BindingDB':
        path = os.getcwd()+'/DBs/BindingDB/tdc_package_preprocessing/BindingDB_max_affinity.tsv'
    elif dataset == 'Davis_et_al':
        path = os.getcwd()+'/DBs/Davis_et_al/tdc_package_preprocessing/DAVIS_et_al.tsv'
    elif dataset == 'BIOSNAP':
        path = os.getcwd()+'/DBs/BIOSNAP/ChG-Miner_miner-chem-gene/ChG-Miner_miner-chem-gene.tsv'
    elif dataset == 'DrugBank':
        path = os.getcwd()+'/DBs/DrugBank/DrugBank_DTIs.tsv'

    #
    def getProts():
        if dataset == 'Yamanashi_et_al_GoldStandard':
            return sorted(np.unique(pd.read_csv(path, sep='\t', header=None, names=['Drug','Protein'])['Protein'].values))
        elif dataset == 'BindingDB' or dataset == 'BIOSNAP' or dataset == 'Davis_et_al' or dataset == 'DrugBank':
            if dataset == 'BindingDB':
                protsids = sorted(np.unique(pd.read_csv(path,sep='\t',index_col=0).Target_ID.values))
            elif dataset == 'BIOSNAP':
                protsids = sorted(np.unique(pd.read_csv(path,delimiter='\t',index_col=0).Gene.values))
            elif dataset == 'DrugBank':
                    protsids = sorted(np.unique(pd.read_csv(path,delimiter='\t',index_col=0).Protein.values))
            elif dataset == 'Davis_et_al':
                protsymbols = sorted(np.unique(pd.read_csv(path,sep='\t',index_col=0).Target_ID.values))
                #FROM GENE SYMBOL TO KEGG/ OR PUBCHEM AND THEN KEGG OR SOMETHING AND EVENTUALLY KEGG.
                BioMART_DB = pd.read_csv(os.getcwd()+'/DBs/BioMART/mart_export_expanded.txt',sep='\t')
                pubchem_genesymbol_to_genename = BioMART_DB[['UniProtKB Gene Name symbol','UniProtKB Gene Name ID']].dropna()
                symbol_gene_dd = dict(zip(pubchem_genesymbol_to_genename['UniProtKB Gene Name symbol'],pubchem_genesymbol_to_genename['UniProtKB Gene Name ID']))
                protsids = sorted(set([symbol_gene_dd[x] for x in protsymbols if x in symbol_gene_dd]))

            print(f'Protein in pubchem -> {len(protsids)}')
            req = r.get('http://rest.kegg.jp/conv/hsa/uniprot')
            print(f'Request status -> {req}')
            pubchem_kegg_conv = [line.decode('utf-8').split('\t') for line in req.iter_lines()]
            pubchem_kegg_conv_dd = {conv_pair[0][3:]:conv_pair[1] for conv_pair in pubchem_kegg_conv}
            #
            prots_in_kegg = sorted([pubchem_kegg_conv_dd[x] for x in protsids if x in pubchem_kegg_conv_dd])
            print(f'Protein in kegg -> {len(prots_in_kegg)}')
            return prots_in_kegg
    #
    def parsehsa(hsastr):
        hsa_l = len(hsastr)
        parse_pos = min([i for i in range(hsa_l) if hsastr[i].isdigit()])
        #
        return hsastr[0:parse_pos] + ':' + hsastr[parse_pos:]
    #
    #Create empty folder
    if not os.path.isdir(path_to_gsymbol):
        os.mkdir(path_to_gsymbol)

    for geneID in tqdm(getProts(),desc='Saving genes symbols'):
        geneID_path = path_to_gsymbol+'/'+geneID+'.txt'
        #if file does not exists
        if not os.path.isfile(geneID_path):
            print(geneID)
            if dataset == 'BindingDB' or dataset == 'Davis_et_al' or dataset == 'BIOSNAP' or dataset == 'DrugBank':
                symbol_v = getSymbol_Kegg(geneID)
            else:
                symbol_v = getSymbol_Kegg(parsehsa(geneID))
            with open(geneID_path, "w") as f:
                    f.writelines('\n'.join(symbol_v))

#getsymbols
dataset, subdataset, path_to_gsymbol = sys.argv[1], sys.argv[2], sys.argv[3]
getSymbol(dataset, subdataset, path_to_gsymbol)