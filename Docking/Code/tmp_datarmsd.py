import pandas as pd
import os

PKLS_PATH = '../Data/pkls'

prot_rmsd = pd.read_pickle(os.path.join(PKLS_PATH, 'RMSD_full_matrix.pkl')) #
prot_annot = pd.read_pickle(os.path.join(PKLS_PATH, 'proteins_annot_crys.pkl'))


new_annot = pd.DataFrame(prot_rmsd.index, columns=['id'])

new_annot['length'] = new_annot.id.map(dict(zip(prot_annot.id, prot_annot.length)))
new_annot['ec'] = new_annot.id.map(dict(zip(prot_annot.id, prot_annot.ec)))
new_annot['Class'] = new_annot.id.map(dict(zip(prot_annot.id, prot_annot.Class)))
new_annot['source'] = new_annot.id.map(dict(zip(prot_annot.id, prot_annot.source)))


#new_annot.to_csv(os.path.join(out_path, f'prots_annot_{key}.csv'), index = False, sep=";")
#print("")


prot_rmsd.to_csv(os.path.join('../Data', f'full_rmsd.csv'), index = False, sep=";")
new_annot.to_csv(os.path.join('../Data', f'full_annotated_proteins.csv'), index = False, sep=";")