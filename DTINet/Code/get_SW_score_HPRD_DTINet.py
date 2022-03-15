import argparse
import logging
import os, sys, uuid
import requests
import pandas as pd
from re import search
from tqdm import tqdm
from itertools import repeat
import subprocess as sp
import multiprocessing as mp
from shutil import rmtree
from sklearn.preprocessing import MinMaxScaler

import numpy as np

wdir = '../Data/DrugBank'

def getamino_uniprot(proteinID):
    r = requests.get(f'https://www.uniprot.org/uniprot/{proteinID}.fasta')
    aminoseq = ''.join(r.text.split('\n')[1:])
    return proteinID, aminoseq

def get_SW_score(pair1, pair2):
	global PATH
	target1, seq1 = pair1
	target2, seq2 = pair2
	fasta1 = os.path.join(PATH, target1.replace(':', '_')+'.fasta')
	if not os.path.exists(fasta1):
		fasta1 = write_fasta(PATH, target1, seq1)
	fasta2 = os.path.join(PATH, target2.replace(':', '_')+'.fasta')
	fasta2 = write_fasta(PATH, target2, seq2)
	result_ID = str(uuid.uuid4())
	result_file = os.path.join(PATH, result_ID+'_results.txt')
	args = ['/home/margaret/data/gserranos/REST_API_embl/EMBOSS-6.6.0/emboss/water', 
			'-asequence', fasta1 , '-bsequence', fasta2, 
			'-gapopen', '10.0', '-gapext', '0.5', 
			'-outfile', result_file]
	try:
		_ = sp.check_call(args, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
		score = extract_score(result_file)
		os.remove(result_file)
		return score
	except:
		print(target1, target2)



def check_and_create_fasta(target, seq):
	global PATH
	fasta1 = os.path.join(PATH, target.replace(':', '_')+'.fasta')
	if not os.path.exists(fasta1):
		fasta1 = write_fasta(PATH, target, seq)
	return fasta1


targets = np.loadtxt(os.path.join(wdir, 'protein.txt'), dtype=str).tolist() 


 
target_seqs = list(map(getamino_uniprot, tqdm(targets)))


all_SmithWaterman = []
for pair1 in tqdm(target_seqs):
    tmp = []
    if not pair1[1]:
        logging.info(f'No sequence for {pair1[0]}')
        continue
    tmp.extend(repeat(pair1, len(target_seqs)))
    with mp.Pool(processes=mp.cpu_count()-5) as pool:
        results = pool.starmap(get_SW_score, zip(tmp, target_seqs))
    all_SmithWaterman.append(results)

targets = [ target for target, _ in target_seqs ]
SmithWaterman_arr = pd.DataFrame(all_SmithWaterman,columns=targets,index=targets)