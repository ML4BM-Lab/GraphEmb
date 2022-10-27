import pandas as pd
import os
import sys
from tqdm import tqdm

#Define dataset path
dataset_path = os.path.join('NeoDTI','Data')

def evaluate_dataset(dir):

    def get_info_from_twoline(files, dataset_matrices):

        #dti will be always just one matrix
        dti = [matrices for f in files for matrices in dataset_matrices if f in matrices][0]

        #read file
        df = pd.read_csv(os.path.join(dir,dti), sep='\t')

        drugs, proteins, edges = len(set(df.iloc[:,0].values)), len(set(df.iloc[:,1].values)), df.shape[0]

        return drugs, proteins, edges, list(map(str,set(df.iloc[:,0].values))), list(map(str,set(df.iloc[:,1].values)))

    def get_info_from_adjmat(files, dataset_matrices, count_zero = False):

        #get adjmats
        adjmats = [matrices for f in files for matrices in dataset_matrices if f in matrices]

        tot_edges = 0
        zero_rows = 0

        for adjmat in tqdm(adjmats, 'Loading adjacency matrices'):

            #read file
            df = pd.read_csv(os.path.join(dir, adjmat), delim_whitespace=True, header = None)
            
            #shape
            tot_rows, tot_columns = df.shape

            if count_zero:
                zero_rows = ((df != 0).sum() == 0).sum() # side effects with all 0s should not be counted as nodes
                tot_edges -= (df == 0).sum().sum()
            
            tot_edges += (df.shape[0] * df.shape[1])

        return tot_rows,  tot_columns - zero_rows, tot_edges

    #Get all matrices
    dataset_matrices = os.listdir(dir)

    ## DTIs
    DTIs_words = ['DTI']
    dti_drug, dti_prots, dti_edges, _, _ = get_info_from_twoline(DTIs_words, dataset_matrices)

    ## SIDE EFFECTS 
    side_effects_words = ['_se.txt']
    _, side_effects_se , side_effects_edges = get_info_from_adjmat(side_effects_words, dataset_matrices, count_zero = True)

    ## DISEASES
    drug_disease_words = ['drug_disease.txt']
    _, drug_diseases , drug_disease_edges = get_info_from_adjmat(drug_disease_words, dataset_matrices)
    prot_disease_words = ['protein_disease.txt']
    _, prot_diseases , prot_disease_edges = get_info_from_adjmat(prot_disease_words, dataset_matrices)
    diseases = max(drug_diseases,prot_diseases)

    ## DRUG PROPERTIES
    drug_properties_words = ['Similarity_Matrix_Drugs']
    _, _ , drug_properties_edges = get_info_from_adjmat(drug_properties_words, dataset_matrices)

    ## PROTEIN PROPERTIES
    protein_properties_words = ['Similarity_Matrix_Proteins']
    _, _ , protein_properties_edges = get_info_from_adjmat(protein_properties_words, dataset_matrices)

    ## DDI
    DDI_words = ['drug_drug']
    _, _ , DDI_edges = get_info_from_adjmat(DDI_words, dataset_matrices)

    ##PPI
    PPI_words =  ['protein_protein']
    _, _ , PPI_edges = get_info_from_adjmat(PPI_words, dataset_matrices)

    """
    Total number of nodes
    """
    
    print(f"Drugs: {dti_drug}, Proteins: {dti_prots}, side effects: {side_effects_se}, diseases: {diseases}, total: {dti_drug+dti_prots+side_effects_se+diseases} ")
    
    """
    Total number of edges
    """
    print(
          f"DTI: {dti_edges}, drug-side_effects: {side_effects_edges}, drug-diseases: {drug_disease_edges}, prot-diseases: {prot_disease_edges}\
          drug-drug: {drug_properties_edges + DDI_edges}, protein-protein: {protein_properties_edges+PPI_edges},\
          total: {dti_edges+side_effects_edges+drug_disease_edges+prot_disease_edges+drug_properties_edges+protein_properties_edges+DDI_edges+PPI_edges}")

for dataset in tqdm(os.listdir(dataset_path), 'Going through each dataset'):

    print(f'Evaluating dataset {dataset}')

    if dataset == 'Yamanashi_et_al_GoldStandard':

        for subdataset in os.listdir(os.path.join(dataset_path, 'Yamanashi_et_al_GoldStandard')):

            print(f'Evaluating subdataset {subdataset}')

            evaluate_dataset(os.path.join(dataset_path, 'Yamanashi_et_al_GoldStandard',subdataset))

    else:

        evaluate_dataset(os.path.join(dataset_path, dataset))
   
