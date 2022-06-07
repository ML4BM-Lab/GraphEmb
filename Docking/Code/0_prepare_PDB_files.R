library(bio3d)

DATA_PATH = '../PDB_FILES'
list_files = list.files(DATA_PATH)

pdb_files <- paste(DATA_PATH, list_files, sep="/")

prot = pdb_files[1]

pdb <- read.pdb(prot)

cpdb <- clean.pdb(pdb)

attributes(pdb)

attributes(cpdb)


length(pdb$atom)

length(cpdb$atom)


atom.select(pdb, chain='-A')