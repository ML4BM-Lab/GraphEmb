import os
import numpy as np
import pandas as pd
from tqdm import tqdm

def read_file(PATH):
	with open(PATH, 'r') as infl:
		lines = infl.readlines()
	return lines


# BIOSNAP
DRUG = 'DB'

PATH = '/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom'
all_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(PATH) for f in filenames ]

matrices = [fl for fl in all_files if not fl.endswith('DTI.txt')]
matrices = [fl for fl in matrices if not fl.endswith('allTsim_files.txt')]
matrices = [fl for fl in matrices if not fl.endswith('allDsim_files.txt')]



all_data = {}
drug_dict = {}
prot_dict = {}
for fl in tqdm(matrices):
	lines = read_file(fl)
	lines = [ln.split() for ln in lines]
	if len(lines[0][0]) > 1000:
		lines = [ln[0].split(',') for ln in lines]
	# add the header of the matrix
	header = lines[0]
	header = [hd for hd in header if hd]
	header = [hd.strip('\"') for hd in header]
	header = [hd.replace('.0', '') for hd in header]
	header = [hd for hd in header if hd != 'data']
	#remove the header
	lines = lines[1:]
	line_id  = [ln[0] for ln in lines]
	line_id = [ln.replace('"', '') for ln in line_id]
	line_id = [ln.replace('.0', '') for ln in line_id]
	data = [ln[1:] for ln in lines]
	# assign to prot
	if not header[0].isnumeric():
		header = [hd.replace(':', '') for hd in header]
		prot_dict[fl] = header
	else:
		drug_dict[fl] = header
	if not line_id[0].isnumeric():
		line_id = [ln.replace(':', '') for ln in line_id]
		prot_dict[fl] = line_id
	else:
		drug_dict[fl] = line_id
	all_data[fl] = (header, line_id, np.array(data,dtype=np.float64))

# select only common ids
all_ids = [set(identifiers) for identifiers in drug_dict.values()]
drugSet = set.intersection(*all_ids)
all_ids = [set(identifiers) for identifiers in prot_dict.values()]
protSet = set.intersection(*all_ids)

for fl in tqdm(matrices):
	header, line_idx, data = all_data[fl]
	if not header[0].isnumeric():
		header_2_keep = [True if hd in protSet  else False for hd in header]
	else:
		header_2_keep = [True if hd in drugSet  else False for hd in header]
	if not line_idx[0].isnumeric():
		lines_2_keep  = [True if hd in protSet  else False for hd in line_idx]
	else:
		lines_2_keep  = [True if hd in drugSet  else False for hd in line_idx]
	column_names =  [ header[i]  for i in range(0,len(header_2_keep))  if header_2_keep[i]]
	row_names  = [ line_idx[i]  for i in range(0,len(lines_2_keep))  if lines_2_keep[i]]
	filtered_data = data[np.ix_(lines_2_keep, header_2_keep)]
	df = pd.DataFrame(filtered_data, columns=column_names, index=row_names)
	df.to_csv(fl , sep='\t', encoding='utf-8', doublequote=False)


with open('/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/DTI.txt', 'r') as infl:
	lines = infl.readlines()
drugprot_set = drugSet.union(protSet)
P_count=[]
D_count=[]
with open('/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom/DTI.txt', 'w') as outfl:
	for ln in lines:
		tmp = ln.split()
		if set(tmp).issubset(drugprot_set):
			P_count.append(tmp[0])
			D_count.append(tmp[1])
			_ = outfl.write('\t'.join(tmp)+'\n')
