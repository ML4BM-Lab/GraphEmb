import enum
import time
from tqdm import tqdm
import logging
import glob
from pymol import cmd
import numpy as np
import pandas as pd
import os
import sys
from tqdm import tqdm
from pymol import cmd
from itertools import product
import multiprocessing as mp

def align_all_to_all(object_list=None, selection='name ca', cutoff=2, cycles=5, debug=0, full_matrix=0, method='super', index_i = 0, index_j =0, out_path='/media/scratch_ssd/tmp/'):
  """
  Copyright (c) 2004 Robert L. Campbell (rlc1@queensu.ca)
  Feel free to do whatever you like with this code.

  Aligns all models in a list to all other models in the list

  usage:
    align_all_to_all [object_list][selection][cutoff=2][cycles=5][debug=0][full_matrix=0][method='align']

        where method can be align, super or cealign

        where selection can be used to specify a particular selection
        (e.g. one chain of each protein),

        cutoff and cycles are options passed to the align or super command.

    By default, all objects are aligned to all others using just the C-alpha atoms.
    You can specify a list of objects to use instead with the object_list option.

    Setting debug=1 prints more information to the terminal or external GUI.

    Setting full_matrix=1 prints out the full symmetric matrix, rather than
    simply the top-half matrix


    Example:
      align_all_to_all object_list=name1 name2 name3 name4, selection=c. a & n. ca, full_matrix=1
  """
  cutoff = int(cutoff)
  full_matrix = int(full_matrix)
  cycles = int(cycles)
  debug=int(debug)

  if not object_list:
    object_list = cmd.get_names()
  else:
    object_list = object_list.replace('[','').replace(']','').replace(',',' ').split()

  rmsd = {}
  rmsd_list = []
#  print object_list
  for i in range(len(object_list)):
    for j in range(i+1,len(object_list)):#, desc=f'loop for it: {i+1} of {len(object_list)}'):
      if method == 'align':
        rms = cmd.align('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection),cutoff=cutoff,cycles=cycles)
      elif method == 'super':
        rms = cmd.super('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection),cutoff=cutoff,cycles=cycles)
      elif method == 'cealign':
        rmsdict = cmd.cealign('%s & %s' % (object_list[i],selection),'%s & %s' % (object_list[j],selection))
        rms = [rmsdict['RMSD'],rmsdict['alignment_length'],1,0,0]
      elif method == 'rms_cur':
#        print'mobile: %s & %s' % (object_list[j],selection),'target: %s & %s' % (object_list[i],selection)
        num_atoms = cmd.select('junkselection','%s & %s' % (object_list[j],selection))
        cmd.delete('junkselection')
        rms = [cmd.rms_cur('%s & %s' % (object_list[j],selection),'%s & %s' % (object_list[i],selection)),num_atoms]

      else:
        print("Method: ",method)
        print("only 'align', 'super', 'cealign' and 'rms_cur' are accepted as methods")
        sys.exit(-1)

      rmsd.setdefault(object_list[i],{})[object_list[j]] = rms[0]
      rmsd_list.append((object_list[j],object_list[i],rms[0],rms[1]))
      if debug and method != 'rms_cur':
        print("Alignment of %s to %s:" % (object_list[j],object_list[i]))
        print("     Initial RMS: %6.3f for %d atoms" % (rms[3],rms[4]))
        print("     Final RMS: %6.3f for %d atoms after %d cycles\n" % (rms[0],rms[1],rms[2]))

#  sort RMSD list
  rmsd_list.sort(key=lambda x: x[2])

# loop over dictionary and print out matrix of final rms values
  if debug:
    for object_name in object_list[:-1]:
      print("%s: %s" % (object_name,str(rmsd[object_name])))

    print("\nSorted from best match to worst:")
    for r in rmsd_list:
      print("%s to %s: %6.3f using %d atoms" % r)
    print("")

  #print("%6s" % " ", end=' ')
  if full_matrix:
  # fill in other half of matrix
    with open(f'{out_path}/test_fullmatrix_{index_i}_{index_j}.txt','w') as file_out:
      for i in range(len(object_list)):
        for j in range(i+1,len(object_list)):
          rmsd.setdefault(object_list[j],{})[object_list[i]] = rmsd[object_list[i]][object_list[j]]
        rmsd[object_list[i]][object_list[i]] = 0
      file_out.write("%6s\t" % " ")
      for i in range(len(rmsd)):
        #print("%6s" % object_list[i], end=' ')
        file_out.write("%6s\t" % object_list[i])
      #print("")
      file_out.write("\n")
      for i in range(len(object_list)):
        #print("%6s" % object_list[i], end=' ')
        file_out.write("%6s\t" % object_list[i])
        for j in range(len(object_list)):
          #print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
          file_out.write("%6.3f\t" % (rmsd[object_list[i]][object_list[j]]))
        #print("")
        file_out.write("\n")
  else:
    with open(f'RMSD_triang_matrix__{index_i}_{index_j}.txt','w') as file_out:
      for i in range(len(rmsd)):
        print("%6s" % object_list[i+1], end=' ')
        file_out.write(f"{object_list[i+1]:6s}")
      #print("")
      file_out.write(f"\n")

      for i in range(len(object_list)):
        #print("%6s" % object_list[i], end=' ')
        file_out.write(f"{object_list[i]:6s}")
        for k in range(i):
          #print("%6s" % " ", end=' ')
          file_out.write(f"{' ':6s}")
        for j in range(i+1,len(object_list)):
          #print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
          file_out.write(f"{rmsd[object_list[i]][object_list[j]]:6.3f}")
        #print("")
        file_out.write(f"\n")


def get_chunk_size(list_objects):
    length = len(list_objects)
    if length > 1000:
        chunk_size = 64
    else:
        chunk_size = 32
    return chunk_size



def row_chunk(chunked_index, i, j):
    list_to_calculate = chunked_index[i] + chunked_index[j]
    # load objects
    for object in list_to_calculate:
            cmd.load(object)
    # calculate and save
    align_all_to_all(object_list=None, selection='name ca', cutoff=2, cycles=5, debug=0, full_matrix=1, method='super', index_i = i, index_j =j, out_path=out_folder)
    # Clean selection; safety check
    cmd.delete('all')
    assert len(cmd.get_object_list()) == 0, 'Previous objects not removed'


##########################################
############### START CODE ###############
logging.basicConfig()
logging.getLogger('').setLevel(logging.INFO)


#### Load Objects
list_pdbs = glob.glob('../Data/Clean_from_PDB/*.pdb')
list_afold = glob.glob('../Data/Clean_from_AFold/*.pdb')
list_objects = list_pdbs + list_afold

# Define and prepare chunks
chunk_size = 52 
ncores = chunk_size
chunked_index = [i.tolist() for i in np.array_split(list_objects, chunk_size)]
# adapt for later create the matrix
list_objects_names = [item.replace('../Data/Clean_from_AFold/', '').replace('../Data/Clean_from_PDB/', '').replace('.pdb', '') for item in list_objects]

logging.info(f'Created {len(chunked_index[0])} chunks with size  {chunk_size}')

# Calculate RMSD
# Superimpose is recommended by Pymol: 
# it is more robust than align for proteins with low sequence similarity
print('Calculating RMSD (method = super)...')

out_folder = '/media/scratch_ssd/tmp/rmsd'
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

logging.info(f'will retrieve {len(chunked_index)**2} files')

###  With multiprocessing
a = time.time()
for i in tqdm(range(len(chunked_index))):
    aa = time.time()
    logging.info(f'Working with chunk {i+1}')
    # Iterating in columns first
    with mp.Pool(processes=ncores) as pool: # change here for less cpu!
        chunked_index_tmp = [chunked_index] * len(chunked_index)
        i_tmp = [i] * len(chunked_index)
        pool.starmap(row_chunk, zip(chunked_index_tmp, i_tmp, range(len(chunked_index)))) #   j should be generated
    print('i finished')
    bb = time.time()
    logging.info(f'Time elapsed for chunk {i+1}: {(bb-aa)/60:.4f} min')

b = time.time()
tm = b-a
logging.info(f'Total Time elapsed: {(tm)/60:.4f} min')



########### FOR READING RESULTS
# later revome because both are chunked index
# this is for creating the full matrix
df = pd.DataFrame(index=list_objects_names, columns=list_objects_names)

for i, j in product(range(len(chunked_index)), range(len(chunked_index))):
    load_data = pd.read_csv(f'{out_folder}/test_fullmatrix_{i}_{j}.txt',  header=0, sep ="\t+", index_col=0,  engine='python')
    chunk_rows = load_data.index.tolist()
    chunk_cols = load_data.columns.tolist()
    df.loc[chunk_rows, chunk_cols] = load_data.to_numpy()

PATH_OUT = '../Results/RMSD_full_matrix.pkl'
logging.info(f'saving files in {PATH_OUT}')
df.to_pickle(PATH_OUT)


logging.debug('Removing temporal files...')
from shutil import rmtree
rmtree(out_folder)
