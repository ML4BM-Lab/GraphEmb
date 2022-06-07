library(bio3d)


# http://thegrantlab.org/bio3d/articles/online/pca_vignette/Bio3D_pca.html


# 

ids = sub('\\.pdb$', '', list.files('afold_files'))[1:20]
folder <- 'afold_files'
pdb_files <- paste(folder, list.files('afold_files'), sep="/")[1:20]

files <- pdbsplit(pdb_files, ids)
pdbs <- pdbaln(files, exefile='./muscle3')

res <- rmsd(pdbs, fit=TRUE)
res[1:10]
write.table(res, file='test.tsv', quote=FALSE, sep='\t')