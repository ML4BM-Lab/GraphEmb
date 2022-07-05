import os
import pandas as pd
import numpy as np
import networkx as nx
import logging

## Extended Datasets
logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)


# list_models = ['DTINet', 'EEG-DTI']
# list_datasets = ['BindingDB',  'BIOSNAP',  'Davis_et_al',  'DrugBank',  'NR', 'E', 'GPCR', 'IC']

# model = 'DTINet'
# dataset = 'NR'

def get_path_folder(model, dataset):
    if dataset in ['NR', 'E', 'GPCR', 'IC']:
        folder_path = f'../../{model}/Data/Yamanashi_et_al_GoldStandard/{dataset}'
    else:
        folder_path = f'../../{model}/Data/{dataset}'
    return folder_path

#


def get_data_DTINet(model='DTINet'):
    list_datasets = ['BindingDB',  'BIOSNAP',  'Davis_et_al',  'DrugBank',  'NR', 'E', 'GPCR', 'IC']
    all_data_model = {}
    for dataset in list_datasets:
        logging.debug(f'Working in {model} with {dataset}')
        # setting folder
        folder_path = get_path_folder(model, dataset)
        #
        # loading data
        dtis = pd.read_csv(os.path.join(folder_path, f'final_dtis_{dataset}.tsv'), sep='\t', header=None)
        ppis = pd.read_csv(os.path.join(folder_path, f'edgelist_PPI.tsv'), sep='\t')
        dds = pd.read_csv(os.path.join(folder_path, f'edgelist_drug_drug.tsv'), sep='\t')
        ddis =  pd.read_csv(os.path.join(folder_path, f'edgelist_drug_disease.tsv'), sep='\t')
        dses  =  pd.read_csv(os.path.join(folder_path, f'edgelist_drug_se.tsv'), sep='\t')
        pdis = pd.read_csv(os.path.join(folder_path, f'edgelist_protein_disease.tsv'), sep='\t')
        # setting all drugs
        drugs = list(set([drug for drug,_ in dtis.values]))
        proteins = list(set([protein for _,protein in dtis.values]))
        # get actual number of diseases, se and etc
        dses = dses[dses.DrugBank_ID.isin(drugs)]
        se = dses.se.unique().tolist()
        # diseases
        ddis = ddis[ddis.DrugBankID.isin(drugs)]
        pdis = pdis[pdis.UniprotID.isin(proteins)]
        diseases = list(set(ddis.DiseaseID.unique().tolist() + pdis.DiseaseID.unique().tolist()))
        # protein protein
        ppis = ppis[ppis.Prot1.isin(proteins)]
        ppis = ppis[ppis.Prot2.isin(proteins)]
        # drug drug
        dds = dds[dds.Drug_ID_A.isin(drugs)]
        dds = dds[dds.Drug_ID_B.isin(drugs)]
        num_drug_nodes = len(drugs)
        num_prot_nodes = len(proteins)
        num_dis_nodes = len(diseases)
        num_se_nodes = len(se)
        total_num_nodes = num_drug_nodes + num_prot_nodes + num_dis_nodes + num_se_nodes
        # load graphs
        Gdti = nx.from_edgelist(dtis.values) # drug-target interactiosn
        Gppi = nx.from_edgelist(ppis.values) # protein-protein interactions
        Gdds = nx.from_edgelist(dds.values) # drug-drug interactions
        Gddis = nx.from_edgelist(ddis.values) # drug-disease
        Gdses = nx.from_edgelist(dses.values) # drug-se
        Gpdis = nx.from_edgelist(pdis.values) # drug disease
        # create full graph
        Gfull = nx.Graph()
        Gfull.add_edges_from(Gdti.edges())
        Gfull.add_edges_from(Gppi.edges())
        Gfull.add_edges_from(Gdds.edges())
        Gfull.add_edges_from(Gddis.edges())
        Gfull.add_edges_from(Gdses.edges())
        Gfull.add_edges_from(Gpdis.edges())
        #
        total_edges = len(Gfull.edges())
        p = num_drug_nodes * num_prot_nodes * num_dis_nodes * num_se_nodes
        sparsity =  total_edges / (p-total_edges)
        connected_componnents = nx.number_connected_components(Gfull)
        list_number_degrees = [degree for _, degree in Gfull.degree]
        list_subgraphs_size = [len(sublist) for sublist in list(nx.connected_components(Gfull))]
        #
        #
        dataset_info = {
            'nodes_drugs': num_drug_nodes,
            'nodes_proteins': num_prot_nodes,
            'nodes_diasease': num_dis_nodes,
            'nodes_se': num_se_nodes,
            'total_nodes': total_num_nodes,
            'dti_edges': len(Gdti.edges),
            'ppi_edges': len(Gppi.edges),
            'drdr_edges': len(Gdds.edges),
            'ddis_edges': len(Gddis.edges),
            'dse_edges': len(Gdses.edges),
            'pdis_edges': len(Gpdis.edges),
            'total_edges':  len(Gfull.edges()),
            'sparsity': sparsity,
            'connected_components': connected_componnents,
            'list_number_degrees': list_number_degrees,
            'list_subgraphs_size': list_subgraphs_size
        }
        all_data_model[dataset] = dataset_info
    return all_data_model


dict_data_dtinet = get_data_DTINet(model='DTINet')

dict_data_eegdti = get_data_DTINet(model='EEG-DTI')




for dataset in list(dict_data_dtinet.keys()):
    print(f'Dataset: {dataset} with sparsity {dict_data_dtinet[dataset]["sparsity"]}')



import json

with open('../Results/statistics_extended_dtinet.json', 'w') as outfile:
    json.dump(dict_data_dtinet, outfile)
  

with open('../Results/statistics_extended_eegdti.json', 'w') as outfile:
    json.dump(dict_data_eegdti, outfile)
  

df_data_dtinet = pd.DataFrame.from_dict(dict_data_dtinet)
df_only_data_dtinet = df_data_dtinet.iloc[:-2,:]
df_only_data_dtinet.to_excel('../Results/statistics_extended_dtinet.xlsx')


df_data_eegdti = pd.DataFrame.from_dict(dict_data_eegdti)

df_only_data_eegdti = df_data_eegdti.iloc[:-2,:]
df_only_data_eegdti.to_excel('../Results/statistics_extended_eeg-dti.xlsx')
