import os
import pandas as pd
import numpy as np
from collections import defaultdict as dd
import random as r
import functools as f
from tqdm import tqdm
from sklearn.model_selection import train_test_split
from random import randint
import Model_Sp_Sd_St_split_Improved as splitter


dti_file = './../../DB/Data/Yamanashi_et_al_GoldStandard/E/interactions/e_admat_dgc_mat_2_line.txt'
DTIs = pd.read_csv(dti_file, sep='\t')
DTIs.columns = ['Drug', 'Protein']

sp_splits = splitter.generate_splits(DTIs, mode= 'Sp', subsampling=True, foldnum=10)