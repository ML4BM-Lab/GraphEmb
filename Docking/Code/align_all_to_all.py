#! /usr/bin/env python3
# Copyright (c) 2004 Robert L. Campbell (rlc1@queensu.ca)

import os
import sys
from tqdm import tqdm
from pymol import cmd

def align_all_to_all(object_list=None,selection='name ca',cutoff=2,cycles=5,debug=0,full_matrix=0,method='align'):
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
    for j in tqdm(range(i+1,len(object_list)), desc=f'loop for it: {i+1} of {len(object_list)}'):
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

  print("%6s" % " ", end=' ')
  if full_matrix:
  # fill in other half of matrix
    with open('test_E.txt','w') as file_out:
      for i in range(len(object_list)):
        for j in range(i+1,len(object_list)):
          rmsd.setdefault(object_list[j],{})[object_list[i]] = rmsd[object_list[i]][object_list[j]]
        rmsd[object_list[i]][object_list[i]] = 0
      file_out.write("%6s\t" % " ")
      for i in range(len(rmsd)):
        print("%6s" % object_list[i], end=' ')
        file_out.write("%6s\t" % object_list[i])
      print("")
      file_out.write("\n")
      for i in range(len(object_list)):
        print("%6s" % object_list[i], end=' ')
        file_out.write("%6s\t" % object_list[i])
        for j in range(len(object_list)):
          print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
          file_out.write("%6.3f\t" % (rmsd[object_list[i]][object_list[j]]))
        print("")
        file_out.write("\n")
  else:
    with open('test.txt','w') as file_out:
      for i in range(len(rmsd)):
        print("%6s" % object_list[i+1], end=' ')
        file_out.write(f"{object_list[i+1]:6s}")
      print("")
      file_out.write(f"\n")

      for i in range(len(object_list)):
        print("%6s" % object_list[i], end=' ')
        file_out.write(f"{object_list[i]:6s}")
        for k in range(i):
          print("%6s" % " ", end=' ')
          file_out.write(f"{' ':6s}")
        for j in range(i+1,len(object_list)):
          print("%6.3f" % (rmsd[object_list[i]][object_list[j]]), end=' ')
          file_out.write(f"{rmsd[object_list[i]][object_list[j]]:6.3f}")
        print("")
        file_out.write(f"\n")

    print(f"outfile:{os.path.join(os.getcwd(),'test.txt')}")

cmd.extend('align_all_to_all',align_all_to_all)
