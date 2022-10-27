import os
import numpy as np

# YAMANISHIS

PATH = '/Drug-Target-Interaction-Prediciton-Method/DTIs_node2vec/Input/Custom'
all_files = [os.path.join(dp, f) for dp, dn, filenames in os.walk(PATH) for f in filenames ]

matrices = [fl for fl in all_files if not fl.endswith('DTI.txt')]
matrices = [fl for fl in matrices if not fl.endswith('allTsim_files.txt')]
matrices = [fl for fl in matrices if not fl.endswith('allDsim_files.txt')]


DRUG = 'D'
PROT = 'hsa'

drug_dict = {}
prot_dict = {}
for fl in matrices:
	with open(fl, 'r') as infl:
		lines = infl.readlines()
	#
	lines = [ln.split() for ln in lines]
	# add the header of the matrix
	header = lines[0]
	header = [hd.strip('\"') for hd in header]
	# add the id of the lines
	lines = lines[1:]
	line_id = [ln[0] for ln in lines]
	line_id = [ln.strip('\"') for ln in line_id]
	if header[0].startswith(PROT):
		header = [hd.replace(':', '') for hd in header]
		prot_dict[fl] = header
	else:
		drug_dict[fl] = header
	if line_id[0].startswith(PROT):
		line_id = [ln.replace(':', '') for ln in line_id]
		prot_dict[fl] = line_id
	else:
		drug_dict[fl] = line_id

# select only common ids
all_ids = [set(identifiers) for identifiers in drug_dict.values()]
set.intersection(*all_ids)

all_ids = [set(identifiers) for identifiers in prot_dict.values()]
set.intersection(*all_ids)

for k, v in drug_dict.items():
	print(k)
	print(v[:3])

for k, v in prot_dict.items():
	print(k)
	print(v[:3])