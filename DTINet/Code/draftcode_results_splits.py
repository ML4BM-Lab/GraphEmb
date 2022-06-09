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
    if db_name == 'Davis':
        split_type = file_info[3]
        subsamp = file_info[4]
    else:
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


nr = pd.read_csv('NR/NR_sub.results', sep="\t")
gpcr = pd.read_csv('GPCR/GPCR_sub.results', sep="\t")
e = pd.read_csv('E/E_sub.results', sep="\t")
ic = pd.read_csv('IC/IC_sub.results', sep="\t")
bindingdb = pd.read_csv('BindingDB/BindingDB_sub.results', sep="\t")
biosnap = pd.read_csv('BIOSNAP/BIOSNAP_sub.results', sep="\t")
drugbank = pd.read_csv('DrugBank/DrugBank_sub.results', sep="\t")
davis = pd.read_csv('Davis_et_al/Davis_sub.results', sep="\t")

plt.clf()

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

axs[0].title.set_text('AUROC ; with subsampling')
axs[1].title.set_text('AUPR ; with subsampling')

dblist = [nr, gpcr, e, bindingdb, biosnap, drugbank]
dblist_name = ['NR', 'GPCR', 'E', 'BindingDB', 'BIOSNAP', 'DrugBank']

for i in range(len(dblist)):
    axs[0].scatter(dblist[i].split_type, dblist[i].AUROC, label=dblist_name[i])
    axs[0].plot(dblist[i].split_type, dblist[i].AUROC)
    axs[0].legend()
    axs[1].scatter(dblist[i].split_type, dblist[i].AUPR, label=dblist_name[i])
    axs[1].plot(dblist[i].split_type, dblist[i].AUPR)
    axs[1].legend()



fig.savefig('test_com.pdf')

import seaborn as sns
sns.color_palette("hls", 8)

plt.clf()

dblist_name = ['NR', 'GPCR', 'E', 'IC','BindingDB', 'BIOSNAP', 'DrugBank', 'Davis']

sub = ['sub', 'nosub']

sam_type = 1

fig, axs = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
if sam_type == 0:
    axs[0].title.set_text('AUROC (with subsampling)')
    axs[1].title.set_text('AUPR  (with subsampling)')
elif sam_type == 1:
    axs[0].title.set_text('AUROC (without subsampling)')
    axs[1].title.set_text('AUPR  (without subsampling)')

for db_name in dblist_name:
    print(db_name)
    if db_name == 'Davis':
        db = pd.read_csv(f'{db_name}_et_al/{db_name}_{sub[sam_type]}.results', sep="\t", index_col=0)
    else:
        db = pd.read_csv(f'{db_name}/{db_name}_{sub[sam_type]}.results', sep="\t", index_col=0)
    db.head()
    axs[0].scatter(db.split_type, db.AUROC, label=db_name)
    axs[0].plot(db.split_type, db.AUROC)
    axs[0].legend(loc=3, prop={'size':8})
    axs[1].scatter(db.split_type, db.AUPR, label=db_name)
    axs[1].plot(db.split_type, db.AUPR)
    axs[1].legend(loc=3, prop={'size':8})

fig.savefig(f'comparison_w_{sub[sam_type]}.png', dpi=400)

