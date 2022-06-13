import os
import pandas as pd


#read df
davis_df = pd.read_csv(os.path.join('.','tdc_package_preprocessing','DAVIS_et_al_w_labels.tsv'),sep='\t')

#get only 2 cols
df2line = davis_df.loc[:,['Drug_ID','Target_ID','Label']]
df2line = df2line[df2line.Label == 1].loc[:,['Drug_ID','Target_ID']]

df2line.to_csv(os.path.join('.','tdc_package_preprocessing','DAVIS_et_al_2line.tsv'),sep='\t',header=False,index=False)


# drugs = pd.read_csv(os.path.join('../../../','DDR/Data/Davis_et_al/Formatted/PreSNF/symmat','Form_symmat_preSNF_Davis_et_al_drug_adjmat_FDA_aers_bit.tsv'),sep='\t',index_col=0)
# drug_names = list(map(str,drugs.index))

# prots = pd.read_csv(os.path.join('../../../','DDR/Data/Davis_et_al/Formatted/PreSNF/symmat','Form_symmat_preSNF_Davis_et_al_prot_SmithWaterman_scores_MinMax.tsv'),sep='\t',index_col=0)
# prots_names = list(map(str,prots.index))

##
# drop_row = []
# for i in range(df2line.shape[0]):
#     if not i%300:
#         print(i)
#     if (str(df2line.iloc[[i],0].values[0]) not in drug_names) or (str(df2line.iloc[[i],1].values[0]) not in prots_names):
#         drop_row.append(i)

# #reduce dataset
# df2linereduced = df2line.drop(drop_row)
# df2linereduced.to_csv(os.path.join('.','tdc_package_preprocessing','DAVIS_et_al_2line.tsv'),sep='\t',header=False,index=False)