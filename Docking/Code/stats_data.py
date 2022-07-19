import os
import pandas as pd
import numpy as np
import logging
import networkx as nx
import json
import helper_functions as hf

logging.basicConfig()
logging.getLogger('').setLevel(logging.DEBUG)

# cols = ['dataset', 'char', 'value']
# pd.DataFrame(columns=cols)



#####

odtis = hf.original_dtis()

# Load all DTIS in a dicionary
dict_dfs  = {
             'DrugBank': odtis.drugbank(),
             'Davis_et_al': odtis.davis(),
             'BindingDB': odtis.bindingdb(),
             'BIOSNAP': odtis.biosnap()
             }

dict_dfs.update(odtis.dict_yamanishi())



all_datasets_info = {}

for key in list(dict_dfs.keys()):
        logging.debug(f'Working in {key}')
        dtis = dict_dfs.get(key)
        records_ = dtis.to_records(index=False)
        G = nx.from_edgelist(records_)
        #
        drugs = [drug for drug,_ in records_]
        proteins = [protein for _,protein in records_]
        #
        drugs = list(set(drugs))
        proteins = list(set(proteins))
        num_drug_nodes = len(drugs)
        num_prot_nodes = len(proteins)
        tot_nodes = num_drug_nodes + num_prot_nodes 
        n_edges = dtis.shape[0]
        #
        adj = nx.adjacency_matrix(G)
        adj_mat = adj.todense()
        #sparsity = round(1-(adj_mat.sum()/adj_mat.shape[0]**2), 6)
        sparsity = n_edges / (num_drug_nodes * num_prot_nodes - n_edges)
        sparsity = round(sparsity, 6)
        connected_componnents = nx.number_connected_components(G)
        list_number_degrees = [degree for _, degree in G.degree]
        list_number_degrees_drugs = [degree for node, degree in G.degree if node in drugs]
        list_number_degrees_proteins = [degree for node, degree in G.degree if node in proteins]
        assert len(list_number_degrees) == len(list_number_degrees_drugs) + len(list_number_degrees_proteins), 'invalid list'
        list_subgraphs_size = [len(sublist) for sublist in list(nx.connected_components(G))]
        #
        dataset_info = {
                'nodes_drugs' : num_drug_nodes,
                'nodes_proteins' : num_prot_nodes,
                'total_nodes' : tot_nodes,
                'total_edges' : n_edges,
                'sparsity_ratio' : sparsity,
                'connected_components': connected_componnents,
                'mean_degree_drugs':  np.mean(list_number_degrees_drugs),
                'mean_degree_proteins': np.mean(list_number_degrees_proteins),
                'list_number_degrees' : list_number_degrees,
                'list_num_degree_drugs': list_number_degrees_drugs,
                'list_num_degree_proteins': list_number_degrees_proteins,
                'list_subgraphs_size' : list_subgraphs_size
        }
        all_datasets_info[key] = dataset_info


df_data = pd.DataFrame.from_dict(all_datasets_info)

df_only_data = df_data.iloc[:-4,:]
df_only_data.to_excel('../Results/statistics_datasets.xlsx')

with open('../Results/statistics_only_dtis.json', 'w') as outfile:
    json.dump(all_datasets_info, outfile)
  

