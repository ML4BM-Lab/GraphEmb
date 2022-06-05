from tkinter import E
import pandas as pd
import numpy as np
import os
import glob
import matplotlib.pyplot as plt


list_files = glob.glob('*.out')
# create db
subsamp_data = []
nosubsamp_data = []

for file in list_files:
#file = list_files[0]
    with open(file) as f:
        lines = f.readlines()
        lines = [line.rstrip() for line in lines]
        line_results = lines[-2]
        line_time = lines[-1]
    time = line_time.split(' ')[3]
    line_results = line_results.replace('Mean: ', '').split(',')
    auroc, aupr = line_results
    auroc = float(auroc.replace('AUROC=', ''))
    aupr =  float(aupr.replace('AUPR=', ''))
    # process info
    name = file.replace('log_DTINet_', '')
    name = name.replace('.out', '')
    file_info = name.split('_')
    db_name = file_info[0]
    split_type = file_info[1]
    subsamp = file_info[2]
    print(db_name, split_type, subsamp, auroc, aupr, time)
    if subsamp == 'subsampl':
        add_row = [split_type, auroc, aupr, time]
        subsamp_data.append(add_row)
    else:
        add_row = [split_type, auroc, aupr, time]
        nosubsamp_data.append(add_row)


nr_sub = pd.DataFrame(subsamp_data,columns=['split_type', 'AUROC', 'AUPR', 'Time'])
nr_nosub = pd.DataFrame(nosubsamp_data,columns=['split_type', 'AUROC', 'AUPR', 'Time'])

def order_df(nr_sub):
    sp_order = {'Sp': 1, 'Sd':2, 'St':3}
    nr_sub['order'] = nr_sub.split_type.map(sp_order)
    nr_sub = nr_sub.sort_values(by='order', ignore_index=True)
    nr_sub = nr_sub.drop(columns = ['order'])
    return nr_sub

nr_sub = order_df(nr_sub)
nr_nosub = order_df(nr_nosub)

nr_sub.to_csv(f'{db_name}_sub.results', sep="\t")
nr_nosub.to_csv(f'{db_name}_nosub.results', sep="\t")

plt.clf()
plt.title(f'Results {db_name} with subsampling')
plt.scatter(nr_sub.split_type, nr_sub.AUROC, label='AUROC')
plt.plot(nr_sub.split_type, nr_sub.AUROC)
plt.scatter(nr_sub.split_type, nr_sub.AUPR, label='AUPR')
plt.plot(nr_sub.split_type, nr_sub.AUPR)
plt.legend()
plt.savefig(f'{db_name}_splits.pdf')


######

#
nr = pd.read_csv('NR/nr_sub.results', sep="\t")
nr = order_df(nr)
gpcr = pd.read_csv('GPCR/GPCR_sub.results', sep="\t")
e = pd.read_csv('E/E_sub.results', sep="\t")
bindingdb = pd.read_csv('BindingDB/BindingDB_sub.results', sep="\t")
biosnap = pd.read_csv('BIOSNAP/BIOSNAP_sub.results', sep="\t")

plt.clf()

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

axs[0].title.set_text('AUROC ; with sumsampling')
axs[1].title.set_text('AUPR ; with sumsampling')

dblist = [nr, gpcr, e, bindingdb, biosnap]
dblist_name = ['nr', 'gpcr', 'e', 'bindingdb', 'biosnap']
for i in range(len(dblist)):
    axs[0].scatter(dblist[i].split_type, dblist[i].AUROC, label=dblist_name[i])
    axs[0].plot(dblist[i].split_type, dblist[i].AUROC)
    axs[0].legend()
    axs[1].scatter(dblist[i].split_type, dblist[i].AUPR, label=dblist_name[i])
    axs[1].plot(dblist[i].split_type, dblist[i].AUPR)
    axs[1].legend()

fig.savefig('test_com.pdf')

