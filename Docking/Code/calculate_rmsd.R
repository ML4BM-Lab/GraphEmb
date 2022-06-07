library(bio3d)
library(parallel)

start_time <- Sys.time()


## doing it now for only pdb files as test
# update this code for pdb + alpha fold
# also give statistics from all data etc and what come from where

folder <- "../Data/Clean_from_PDB"
ids <- sub("\\.pdb$", "", list.files(folder))
pdb_files <- paste(folder, list.files(folder), sep = "/")

n_cores <- 1 # change here for both
# check this
print("splitting...")
files <- pdbsplit(pdb_files, ids)
print("aligment...")
pdbs <- pdbaln(files, exefile = "./muscle3")
# alt records selected automatically to A
# can pass ncore = 2 to paralellice
print("rmsd...")
res <- rmsd(pdbs, fit = TRUE)
res[1:10]

end_time <- Sys.time()

write.table(res, file = "all_PDB_prot_rmsd.tsv", quote = FALSE, sep = "\t")


print(end_time - start_time)