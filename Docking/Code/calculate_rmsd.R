library(bio3d)
library(parallel)

start_time <- Sys.time()

folder_pdb <- "../Data/Clean_from_PDB"
ids <- sub("\\.pdb$", "", list.files(folder_pdb))
pdb_files <- paste(folder_pdb, list.files(folder_pdb), sep = "/")

folder_afold <- "../Data/Clean_from_AFold"
ids_a <- sub("\\.pdb$", "", list.files(folder_afold))
pdbs_a <- paste(folder_afold, list.files(folder_afold), sep = "/")

# append afold data; now we have the full list
ids <- append(ids, ids_a)
pdb_files <- append(pdb_files, pdbs_a)

# Number of cores to use
n_cores <- 5
# check this
print("Splitting...")
# overwrite FALSE: PDBs not be read and written if split files already exist
files <- pdbsplit(pdb_files, ids, overwrite = FALSE)
print("Aligment...")
pdbs <- pdbaln(files, ncore = n_cores, exefile = "./muscle3")
#
saveRDS(pdbs, "../Data/pdbs_aln.rds")
# alt records selected automatically to A
print("RMSD calculation...")
res <- rmsd(pdbs, fit = TRUE, ncore = n_cores)

end_time <- Sys.time()

print("Saving results...")
write.table(res, file = "../Data/matrix_rmsd.tsv", quote = FALSE, sep = "\t")

cat("# cores:", n_cores, "\n")
print(end_time - start_time)