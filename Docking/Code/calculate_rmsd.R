library(bio3d)


start_time <- Sys.time()

## doing it now for only pdb files as test
# update this code for pdb + alpha fold
folder <- '../Data/Clean_from_PDB'
ids = sub('\\.pdb$', '', list.files(folder))
pdb_files <- paste(folder, list.files(folder), sep="/")

# check this
files <- pdbsplit(pdb_files, ids)
pdbs <- pdbaln(files, exefile='./muscle3')
# alt records selected automatically to A

res <- rmsd(pdbs, fit=TRUE)
res[1:10]

end_time <- Sys.time()

write.table(res, file='all_PDB_prot_rmsd.tsv', quote=FALSE, sep='\t')


print(end_time - start_time)