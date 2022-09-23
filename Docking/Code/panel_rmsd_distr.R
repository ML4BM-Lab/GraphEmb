
if(!require(ComplexHeatmap)) biocLite("ComplexHeatmap")
#if(!require(circlize)) install.packages('circlize')

library(devtools)
library(circlize)
# better rasterization
library(magick)
library(ggplot2)
library(RColorBrewer)

##################
##################
# LOAD FILES
database_name <- 'NR'

get_data <- function(database_name){
  #
  print(database_name)
  folder_path <- '../Results/data_per_dataset'
  folder_db_path <-  file.path(folder_path, database_name)
  #
  file_prot_sim <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))
  #
  # load
  mat_rmsd_proteins <- read.table(file_prot_sim, sep = ";", header=T)
  mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
  dim(mat_rmsd_proteins)
  #
  ##############
  ## PLOT FIRST HISTOGRAM
  hits_prot <- mat_rmsd_proteins[upper.tri(mat_rmsd_proteins, diag = FALSE)]
  #data_hist <- data.frame(hits_prot)
  #
  # test
  scaler <- function(x){(x-min(x))/(max(x)-min(x))}
  hits_prot <- 1-scaler(hits_prot)
  #data_hist <- data.frame(hits_prot)
return(hist_prot, data_hist)}

#file_fig_hist <- paste0(folder_db_path, '/hist_prot_', database_name, '.pdf')

file_fig_hist <- 'test_panel.pdf'

print(max(mat_rmsd_proteins))

pdf(file = file_fig_hist)
ggplot(data_hist, aes(hits_prot)) +
  ggtitle(database_name)+
  labs(x = "RMSD (Å)") +             
  geom_histogram( bins=60, colour="#0c0f0a", fill="#f18f01")+
  #geom_density()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))   
dev.off()




