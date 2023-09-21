import pandas as pd
import os
import sys
from tqdm import tqdm

#Define dataset path
dataset_path = os.path.join('DDR','Data')

def evaluate_dataset(dataset, dir):

    def get_info_from_twoline(files, dataset_matrices):

        #dti will be always just one matrix
        dti = [matrices for f in files for matrices in dataset_matrices if f in matrices][0]

        #read file
        df = pd.read_csv(os.path.join(dir,dti),sep='\t', header=None)

        drugs, proteins, edges = len(set(df.iloc[:,0].values)), len(set(df.iloc[:,1].values)), df.shape[0]

        return drugs, proteins, edges, list(map(str,set(df.iloc[:,0].values))), list(map(str,set(df.iloc[:,1].values)))

    def get_info_from_adjmat(files, dataset_matrices, reference, count_zero = False):

        #get adjmats
        adjmats = [matrices for f in files for matrices in dataset_matrices if f in matrices]

        tot_rows = []
        tot_columns = []
        tot_edges = 0
        zero_rows = 0

        for adjmat in tqdm(adjmats, 'Loading adjacency matrices'):

            #read file
            df = pd.read_csv(os.path.join(dir, adjmat), delim_whitespace=True, index_col = 0)
            df.index = list(map(str,df.index))
            reference = list(set(df.index).intersection(set(reference)))

            #filter by row
            df_rowfilter = df.loc[reference,:]

            #check if columns apply
            if df.columns[0] in reference:
                df_columnfilter = df_rowfilter.loc[:,reference]
            else:
                df_columnfilter = df_rowfilter.copy()

            tot_rows += df_columnfilter.index.tolist()
            tot_columns += df_columnfilter.columns.tolist()

            if count_zero:
                zero_rows = ((df_columnfilter != 0).sum() == 0).sum() # side effects with all 0s should not be counted as nodes
                tot_edges -= (df_columnfilter == 0).sum().sum()
            
            tot_edges += (df_columnfilter.shape[0] * df_columnfilter.shape[1])

        return len(set(tot_rows)),  len(set(tot_columns)) - zero_rows, tot_edges

    #Get all matrices
    dataset_matrices = os.listdir(dir)

    ## DTIs
    DTIs_words = ['twoline', '2line', '2_line']
    dti_drug, dti_prots, dti_edges, dti_drug_names, dti_prots_names = get_info_from_twoline(DTIs_words, dataset_matrices)

    if dataset in ['NR','IC','E','GPCR']:
        dti_prots_names = [x[0:3]+':'+x[3:] for x in dti_prots_names]

    ## SIDE EFFECTS 
    side_effects_words = ['_drug_FDA_aers_bit.tsv', '_drug_FDA_aers_freq.tsv']
    side_effects_drugs, side_effects_se , side_effects_edges = get_info_from_adjmat(side_effects_words, dataset_matrices, dti_drug_names, count_zero = True)
    #only edges
    side_effect_only_edges_words = ['_drug_SIDER_SideEffect.tsv']
    _, _ , side_effects_only_edges = get_info_from_adjmat(side_effect_only_edges_words, dataset_matrices, dti_drug_names)

    ## DRUG PROPERTIES
    drug_properties_words = ['_drug_Rchemcpp']
    drug_properties_drugs, _ , drug_properties_edges = get_info_from_adjmat(drug_properties_words, dataset_matrices, dti_drug_names)

    ## PROTEIN PROPERTIES
    protein_properties_words = ['_prot_mismatch', '_prot_SmithWaterman', '_prot_spectrum', '_prot_BioGrid']
    protein_properties_prots, _ , protein_properties_edges = get_info_from_adjmat(protein_properties_words, dataset_matrices, dti_prots_names)

    #GO
    go_words = ['GO']
    go_proteins, _ , go_edges = get_info_from_adjmat(go_words, dataset_matrices, dti_prots_names)

    print(f'\nDTIs -> drugs: {dti_drug}, prots: {dti_prots}, edges: {dti_edges}')
    print(f'Side Effects -> drugs: {side_effects_drugs}, side effects: {side_effects_se}, edges: {side_effects_edges}')
    print(f'Side effects, only edges -> edges: {side_effects_only_edges}')
    print(f'Drug Properties -> drugs: {drug_properties_drugs}, edges: {drug_properties_edges}')
    print(f'Protein Properties -> proteins: {protein_properties_prots}, edges: {protein_properties_edges}')
    print(f'GO -> proteins: {go_proteins}, edges: {go_edges}\n')

    # ----------- #
    print(f'\nDTIs -> drugs: {dti_drug}, prots: {dti_prots}, edges: {dti_edges}')
    print(f'Drug-Drug edges: {drug_properties_edges + side_effects_only_edges}')
    print(f'Protein-Protein edges: {protein_properties_edges + go_edges}')

for dataset in tqdm(os.listdir(dataset_path), 'Going through each dataset'):

    print(f'Evaluating dataset {dataset}')

    if dataset == 'Yamanashi_et_al_GoldStandard':

        for subdataset in os.listdir(os.path.join(dataset_path, 'Yamanashi_et_al_GoldStandard')):

            print(f'Evaluating subdataset {subdataset}')

            evaluate_dataset(subdataset, os.path.join(dataset_path, 'Yamanashi_et_al_GoldStandard',subdataset))

    else:

        evaluate_dataset(dataset, os.path.join(dataset_path, dataset))
   
