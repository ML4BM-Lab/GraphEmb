
if(!require(ComplexHeatmap)) biocLite("ComplexHeatmap")
#if(!require(circlize)) install.packages('circlize')

library(devtools)
library(circlize)
# better rasterization
library(magick)
library(ggplot2)
library(RColorBrewer)



# Proteins


get_matrix_proteins <- function(database_name){
    print(database_name)
    folder_path <- '../Results/data_per_dataset'
    folder_db_path <-  file.path(folder_path, database_name)
    file_prot_sim <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))
    mat_rmsd_proteins <- read.table(file_prot_sim, sep = ";", header=T)
    mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
    mat_rmsd_proteins
    }

get_annot_proteins <- function(database_name){
    folder_path <- '../Results/data_per_dataset'
    folder_db_path <-  file.path(folder_path, database_name)
    file_prot_annot <- file.path(folder_path, database_name, paste0('prots_annot_', database_name, '.csv'))
    annots_prot <- read.table(file_prot_annot, header = TRUE, sep = ";", quote = "")
}


database_name = 'NR'

mat_rmsd_proteins <- get_matrix_proteins(database_name)
annots_prot <- get_annot_proteins(database_name)
