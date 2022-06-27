import pandas as pd
import numpy as np
import os
import sys
import functools as ft
from tqdm import tqdm
from numpy import savetxt

def visualize(dataset,subdataset):

    fpath = os.path.join(os.getcwd(),'DDR/Data',dataset,subdataset)

    prot_f = []
    drug_f = []

    for f in os.listdir(fpath):

            if 'prot' in f:
                prot_f.append(f)
            elif 'drug' in f:
                if 'aers' in f and 'adjmat' not in f:
                    pass
                else:
                    drug_f.append(f)

    for f in prot_f:
        
        mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)

        try:
            matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
        except:
            matnames = mat.columns.to_list()

        print(f'{matnames[0:5]} for file {f} (targets)')

    for f in drug_f:

            if 'aers' in f and dataset == 'BindingDB':
                mat = pd.read_csv(os.path.join(fpath,f),sep=',',index_col=0)
            else:
                if 'Rchem' not in f:
                    mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
                else:
                    mat = pd.read_csv(os.path.join(fpath,f),sep=' ',index_col=0)

            try:
                matnames = [drug.replace(':','') for drug in mat.index.to_list()]
            except:
                matnames = mat.index.to_list()

            print(f'{matnames[0:5]} for file {f} (drugs)')

def format_matrices(dataset,subdataset,symmat=True,twolines=True,symmat_nonames=True):

    #define path
    fpath = os.path.join(os.getcwd(),'DDR/Data',dataset,subdataset)
    
    #print(f"path {fpath}")

    def retrieve_matrices():
        
        prot_f = []
        drug_f = []

        #categorize in drug/proteins files
        for f in os.listdir(fpath):
            if 'prot' in f:
                prot_f.append(f)
            elif 'drug' in f:
                if 'aers' in f and 'adjmat' not in f:
                    pass
                else:
                    drug_f.append(f)

        #return the sorted array
        return sorted(prot_f), sorted(drug_f)

    def retrieve_dti():

        dtipath = [f for f in os.listdir(fpath) if '2line' in f][0]

        #read the f
        dtif = pd.read_csv(os.path.join(fpath,dtipath),sep='\t',header=None)
        dtif.columns = ['Drugs','Targets']

        return dtif.Drugs.tolist(), dtif.Targets.tolist()

    def filter_matrices_proteins_symmat(mat_l,savepath,dtil=None):

        names = []

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)

            try:
                matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
            except:
                matnames = mat.columns.to_list()

            print(f'{matnames[0:5]} for file {f} (targets)')
            names.append(matnames)

        #now select the intersection
        inames = sorted(ft.reduce(lambda a,b: set(a).intersection(set(b)),names))

        if dtil:
            inames = sorted(set(inames).intersection(set(dtil)))

        for f in mat_l:
            mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)

            try:
                matnames = [gene.replace(':','') for gene in mat.columns.to_list()]
            except:
                matnames = mat.columns.to_list()
            
            mat.columns = matnames
            mat.index = matnames
            #generate the final mat and save
            fmat = mat.loc[inames,inames]
            print(f"shape of formatted target matrix-> {fmat.shape}")
            fmat.to_csv(os.path.join(savepath,'Form_symmat_preSNF_'+f),sep='\t')

    def filter_matrices_drugs_symmat(mat_l,savepath,dtil=None):

        names = []

        for f in mat_l:

            if 'aers' in f and dataset == 'BindingDB':
                mat = pd.read_csv(os.path.join(fpath,f),sep=',',index_col=0)
            else:
                if 'Rchem' not in f:
                    mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
                else:
                    mat = pd.read_csv(os.path.join(fpath,f),sep=' ',index_col=0)

            try:
                matnames = [drug.replace(':','') for drug in mat.index.to_list()]
            except:
                matnames = mat.index.to_list()

            print(f'{matnames[0:5]} for file {f} (drugs)')
            names.append(matnames)

        #now select the intersection
        inames = sorted(ft.reduce(lambda a,b: set(a).intersection(set(b)),names))

        if dtil:
            print(f"inames type {type(inames[0])}, dtil type {type(dtil[0])}")
            inames = sorted(set(inames).intersection(set(dtil)))

        for f in mat_l:

            if 'aers' in f and dataset == 'BindingDB':
                mat = pd.read_csv(os.path.join(fpath,f),sep=',',index_col=0)
            else:
                if 'Rchem' not in f:
                    mat = pd.read_csv(os.path.join(fpath,f),sep='\t',index_col=0)
                else:
                    mat = pd.read_csv(os.path.join(fpath,f),sep=' ',index_col=0)

            try:
                matnames = [drug.replace(':','') for drug in mat.index.to_list()]
            except:
                matnames = mat.index.to_list()

            mat.index = matnames
            mat.columns = matnames
            
            #intersect names
            fmat = mat.loc[inames,inames]

            print(f"shape of formatted drug matrix-> {fmat.shape}")

            fmat.to_csv(os.path.join(savepath,'Form_symmat_preSNF_'+f),sep='\t')

    def symmat_to_2lines(dirpath,savepath):

        fnames = os.listdir(dirpath)

        for f in tqdm(fnames,desc='Converting files...'):

            fsymmat = pd.read_csv(os.path.join(dirpath,f),sep='\t',index_col=0)
            fsymmatnp = fsymmat.to_numpy()
            
            #getting names and path to save f
            enames = list(map(str,fsymmat.index))
            #enames = fsymmat.index
            fsave = os.path.join(savepath,f)
            fsave = fsave.replace('symmat','2line')

            with open(fsave,'w') as appendf:

                for i in tqdm(range(fsymmat.shape[0]),desc=f'Saving file {f}'):
                    for j in range(i,fsymmat.shape[0]):
                        row = '\t'.join([enames[i],enames[j],str(fsymmatnp[i,j])]) + '\n'
                        appendf.write(row)

    def symmat_nonames(dirpath,savepath):

        fnames = os.listdir(dirpath)

        for f in tqdm(fnames,desc='Converting files...'):

            fsymmat = pd.read_csv(os.path.join(dirpath,f),sep='\t',index_col=0)
            fsymmatnp = np.round(fsymmat.to_numpy(),2)
            savetxt(os.path.join(savepath,f[:-4]+'_nonames.csv'),fsymmatnp,delimiter=',',fmt='%1.1f')

        #finally, generate a list containing all these names
        fnames = os.listdir(savepath)

        with open(os.path.join(savepath,'snflist_drugs.txt'),'w') as snf:
            snf.write('\n'.join([f for f in fnames if 'drug' in f]))

        with open(os.path.join(savepath,'snflist_prots.txt'),'w') as snf:
            snf.write('\n'.join([f for f in fnames if 'prot' in f]))

    def check_SNF(dirpath):
        fnames = os.listdir(dirpath)

        from sklearn.metrics.pairwise import euclidean_distances
        import sklearn

        arr1 = np.loadtxt(os.path.join(dirpath,fnames[0]), delimiter=',')
        arr2 = np.loadtxt(os.path.join(dirpath,fnames[1]), delimiter=',')

        sklearn.metrics.pairwise.pairwise_distances(arr1,arr2,metric='euclidean')

    def post_SNF(dirpath, savepath):

        #drugs
        drugpath = [os.path.join(savepath,f) for f in os.listdir(savepath) if 'drug' in f][0]

        drugfiles = pd.read_csv(drugpath,header=None).values.flatten().tolist()

        for drugfile in drugfiles:
            pass

        #prots
        protpath = [os.path.join(savepath,f) for f in os.listdir(savepath) if 'prot' in f][0]

        with open(protpath,'r') as protf:
            protf.readline()

    prot_f, drug_f = retrieve_matrices()
    #dtidrug, dtiprot = retrieve_dti()

    #print(f'Matrices for proteins {prot_f}')
    #print(f'Matrices for drugs {drug_f}')

    folder_path = os.path.join(fpath,'Formatted/')
    #print(f"folder path {folder_path}")

    if not os.path.isdir(folder_path):
        #create the folder
        os.mkdir(os.path.join(fpath,'Formatted'))

    #---------------------- PreSNF --------------------- #

    folder_path_preSNF = os.path.join(folder_path,'PreSNF/')
    if not os.path.isdir(folder_path_preSNF):
        os.mkdir(folder_path_preSNF)

    folder_path_preSNF_2line = os.path.join(folder_path,'PreSNF/2line/')
    if not os.path.isdir(folder_path_preSNF_2line):
        os.mkdir(folder_path_preSNF_2line)

    folder_path_preSNF_symmat = os.path.join(folder_path,'PreSNF/symmat/')
    if not os.path.isdir(folder_path_preSNF_symmat):
        os.mkdir(folder_path_preSNF_symmat)

    folder_path_preSNF_symmat_nonames = os.path.join(folder_path,'PreSNF/symmat_nonames/')
    if not os.path.isdir(folder_path_preSNF_symmat_nonames):
        os.mkdir(folder_path_preSNF_symmat_nonames)

    #------------------------ PostSNF ---------------------- #

    folder_path_postSNF = os.path.join(folder_path,'PostSNF/')
    if not os.path.isdir(folder_path_postSNF):
        os.mkdir(folder_path_postSNF)

    #------------------- Format matrices ------------------ #

    if symmat:
        print("\ngenerating symmat")
        filter_matrices_proteins_symmat(prot_f,savepath = folder_path_preSNF_symmat)
        filter_matrices_drugs_symmat(drug_f,savepath = folder_path_preSNF_symmat)
    
    if twolines:
        print("\nconverting symmat to twolines")
        symmat_to_2lines(dirpath = folder_path_preSNF_symmat, savepath = folder_path_preSNF_2line)

    if symmat_nonames:
        print("\ngenerating symmat but with no names")
        symmat_nonames(dirpath = folder_path_preSNF_symmat, savepath = folder_path_preSNF_symmat_nonames)

def reformatPost(dataset,subdataset):

    def addl(v):
        for i in range(len(v)):
            if isinstance(v[i],int):
                v[i] = 'A'+str(v[i])

        return v

    #define path
    fpath = os.path.join(os.getcwd(),'DDR/Data',dataset,subdataset)

    #folder paths
    folder_path = os.path.join(fpath,'Formatted/')
    folder_path_postSNF = os.path.join(folder_path,'PostSNF/')

    #get matrices
    postfiles = [f for f in os.listdir(folder_path_postSNF) if '.tsv' in f]

    for f in postfiles:

        #read mat
        mat = pd.read_csv(os.path.join(folder_path_postSNF,f),sep='\t',header=None)
        
        #get data by cols
        zerocol = addl(mat.iloc[:,0].values.tolist())
        onecol = addl(mat.iloc[:,1].values.tolist())

        if mat.shape[1] > 2:
            value = mat.iloc[:,2].values.tolist()
            #generate new dataframe
            newdf = pd.DataFrame([zerocol,onecol,value]).T
        else:
            newdf = pd.DataFrame([zerocol,onecol]).T
        newdf.to_csv(os.path.join(folder_path_postSNF,f),sep='\t',header=False,index=False)

def checkPost(dataset,subdataset,modify=False):

    #define path
    fpath = os.path.join(os.getcwd(),'DDR/Data',dataset,subdataset)

    #folder paths
    folder_path = os.path.join(fpath,'Formatted/')
    folder_path_postSNF = os.path.join(folder_path,'PostSNF/')

    #get matrices
    postfiles = [f for f in os.listdir(folder_path_postSNF) if 'Form' in f]
    dtimat = [f for f in os.listdir(folder_path_postSNF) if 'Form' not in f and f[0:1] == 'B']

    for f in dtimat:

        mat = pd.read_csv(os.path.join(folder_path_postSNF,f),sep='\t',header=None)

        drugs = set(mat.iloc[:,0].values.tolist())
        prots = set(mat.iloc[:,1].values.tolist())

        len_drugs = len(drugs)
        len_prots = len(prots)

        len_adjmat_drugs = len_drugs * len_drugs / 2
        len_adjmat_prots = len_prots * len_prots / 2

        print(f"Shape of dti mat {mat.shape}, shape of drugs {len_drugs} ({len_adjmat_drugs} adjmat), shape of proteins {len_prots} ({len_adjmat_prots} adjmat) ")

    def computecount(matrix,totallen,el,name='drugs'):

        elem1 = set(matrix.iloc[:,0].values.tolist())
        elem2 = set(matrix.iloc[:,1].values.tolist())
        elemt = elem1.union(elem2)

        print(f"adjmat {name} {len(elemt)} against total {name} {totallen}, with an intersection of {len(elemt.intersection(el))}")

        return elemt

    for f in postfiles:
        
        mat = pd.read_csv(os.path.join(folder_path_postSNF,f),sep='\t',header=None)
        print(f"Shape of adjmat {mat.shape}")

        if 'drug' in f:
            drug_names = computecount(mat,len_drugs,drugs,'drugs')
        elif 'prot' in f:
            prots_names = computecount(mat,len_prots,prots,'prots')

    if modify:

        #compute intersect of drugs and genes
        drug_intersection = drugs.intersection(drug_names)
        prot_intersection = prots.intersection(prots_names)

        for f in dtimat:
            drop_row = []
            mat = pd.read_csv(os.path.join(folder_path_postSNF,f),sep='\t',header=None)
            for i in range(mat.shape[0]):
                if not i%300:
                    print(i)
                if (mat.iloc[[i],0].values[0] not in drug_intersection) or (mat.iloc[[i],1].values[0] not in prot_intersection):
                    drop_row.append(i)

            #reduce dataset
            df2linereduced = mat.drop(drop_row)
            #print(f"new shape {df2linereduced.shape}")
            df2linereduced.to_csv(os.path.join(folder_path_postSNF,f),sep='\t',header=False,index=False)   

        checkPost(dataset,subdataset,modify=False)

#get sys
try:
    arg1 = sys.argv[1]
    try: 
        arg2 = sys.argv[2]
    except:
        arg2 = ''
except:
    arg1, arg2 = '','' 

dataset, subdataset = arg1, arg2

#visualize(dataset,subdataset)
format_matrices(dataset,subdataset)
#reformatPost(dataset,subdataset)
#checkPost(dataset,subdataset)