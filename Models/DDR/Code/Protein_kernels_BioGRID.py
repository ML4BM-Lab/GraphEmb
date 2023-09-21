##call python from guille's path
##/home/sevastopol/anaconda3/bin/python

import pandas as pd
import os
import numpy as np
from tqdm import tqdm
import requests as r
import sys
import math as m
from igraph import *

def generate_protein_kernel_BioGRID(savepath,protpath):

    def cmp_gsymbols(cmp_dict,f1,f2):
        #
        #A->B 
        a_to_b_v1 = cmp_dict[f1[:-4]]['v1']
        a_to_b_v2 = cmp_dict[f2[:-4]]['v2']
        #
        #B->A, now f1 is on 
        b_to_a_v1 = cmp_dict[f2[:-4]]['v1']
        b_to_a_v2 = cmp_dict[f1[:-4]]['v2']
        #
        if set(a_to_b_v1).intersection(set(a_to_b_v2)) or set(b_to_a_v1).intersection(set(b_to_a_v2)):
            return True
        else:
            return False

    def generate_cmp_dict(dfiles,rf1,rf2):
        #
        def gen_dict(var):
            vardd = {}
            for i,x in enumerate(var):
                if x in vardd:
                    vardd[x].append(i)
                else:
                    vardd[x] = [i]
            return vardd
        #
        rf1dd = gen_dict(rf1)
        rf2dd = gen_dict(rf2)
        #init cmp dict
        cmp_dict = {}
        #
        #fill it using files
        for f in tqdm(dfiles):
            #
            #print(f)
            try:
                symbol_l = pd.read_csv(protpath+'/'+f,header=None).values.flatten().tolist()
            except:
                cmp_dict[f[:-4]] = {'v1':[],'v2':[]}
                continue
            #
            a_to_b_v1 = []
            a_to_b_v2 = []
            #
            if len(symbol_l):
                #
                a_to_b_v1 = [rf1dd[x] for x in symbol_l if x in rf1dd]
                a_to_b_v2 = [rf2dd[x] for x in symbol_l if x in rf2dd]
                #
                if len(a_to_b_v1):
                    a_to_b_v1 = np.hstack(a_to_b_v1)
                #
                if len(a_to_b_v2):
                    a_to_b_v2 = np.hstack(a_to_b_v2)
            #   
            cmp_dict[f[:-4]] = {'v1':a_to_b_v1,'v2':a_to_b_v2}
        #
        return cmp_dict

    print('Generating symbol files')

    #READ THE WHOLE DATASET
    BioGridDB = pd.read_csv(os.getcwd()+'/DBs/BioGRID/HeaderLess_BIOGRID-ALL-4.4.207.tab.txt',sep='\t')

    #TAKE THE OFFICIAL SYMBOL AND LOOK IF WE HAVE IT
    dfiles = sorted(os.listdir(protpath))
    newdfiles = []

    #check empties
    for f in tqdm(dfiles):
        try:
            _ = pd.read_csv(protpath+'/'+f,header=None).values.flatten().tolist()
            newdfiles.append(f)
        except:
            pass
    #
    print('from ',len(dfiles),' to ',len(newdfiles))
    newdfiles = dfiles

    #build the matrix
    PPI_mat = np.zeros(shape=(len(dfiles),len(dfiles)))

    ##
    official_symbolA = np.array(BioGridDB['OFFICIAL_SYMBOL_A'].values.tolist())
    official_symbolB = np.array(BioGridDB['OFFICIAL_SYMBOL_B'].values.tolist())

    print('Comparing symbols from dataset with symbols in BioGrid DB')
    ##dict
    gsymbol_cmp_dict = generate_cmp_dict(dfiles,official_symbolA,official_symbolB)

    print('Building the adjacency matrix')
    #only triangular
    for i1,f1 in enumerate(tqdm(dfiles)):
        for i2,f2 in enumerate(dfiles[i1:]):
            #select rows 
            PPI_mat[i1,i1+i2] = cmp_gsymbols(gsymbol_cmp_dict,f1,f2)
            PPI_mat[i1+i2,i1] = PPI_mat[i1,i1+i2]
                
    print('Building the graph from adjacency matrix')
    #build the graph from the adjacency matrix
    g = Graph(directed=False).Adjacency(PPI_mat.tolist())

    PPI_mat_shd = np.zeros(shape=(len(dfiles),len(dfiles)))

    print('Computing shortest hop distance')
    for i1,f1 in enumerate(tqdm(dfiles)):
        for i2,f2 in enumerate(dfiles[i1:]):
            #compute shortest path
            hop_shortest_distance = g.shortest_paths(i1,i1+i2)[0][0]
            # if hop_shortest_distance:
            #     print(hop_shortest_distance)
            #use formula
            hop_shortest_distance = 0.9 * m.exp(-1*hop_shortest_distance)
            #S(p,p') = A * e^b*D(p,p')
            #A=0.9, b=1 and D(p,p') is the shortest hop distance
            #assign 
            PPI_mat_shd[i1,i1+i2] = hop_shortest_distance
            PPI_mat_shd[i1+i2,i1] = hop_shortest_distance

    PPI_mat_shd_pd = pd.DataFrame(PPI_mat_shd,index=[f[:-4] for f in dfiles],columns=[f[:-4] for f in dfiles])
    print(f'Saving matrix, with shape {PPI_mat_shd_pd.shape}')
    PPI_mat_shd_pd.to_csv(savepath+'BioGrid_SHD.tsv',sep='\t')
    
#generate protein kernel using BioGRID DB and shortest hop distance metric
generate_protein_kernel_BioGRID(sys.argv[1],sys.argv[2])