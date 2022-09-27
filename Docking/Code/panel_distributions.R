
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

scaler <- function(x){(x-min(x))/(max(x)-min(x))}

get_data_rmsd <- function(database_name){
  print(database_name)
  folder_path <- '../Results/data_per_dataset'
  folder_db_path <-  file.path(folder_path, database_name)
  file_prot_sim <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))
  # load
  print(file_prot_sim)
  mat_rmsd_proteins <- read.table(file_prot_sim, sep = ";", header=T)
  mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
  print( dim(mat_rmsd_proteins))
  hits_prot <- mat_rmsd_proteins[upper.tri(mat_rmsd_proteins, diag = FALSE)]
  hits_prot <- scaler(hits_prot)
  tit <- paste0(database_name, ' max.val. ', max(mat_rmsd_proteins) )
  data_hist <- data.frame(hits_prot)
    ggplot(data_hist, aes(hits_prot)) +
    ggtitle(tit)+
    labs(x = "RMSD (Å)") +             
    geom_histogram(bins=60, colour="#0c0f0a", fill="#E03B29")+
    #geom_density()+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size=8), 
          axis.title = element_text(size=5),
          axis.text = element_text(size=5))   
}

#data_hist <- data.frame(hits_prot)
#
list_datasets = c('NR', 'IC', 'GPCR', 'E', 'DrugBank', 'Davis_et_al', 'BindingDB', 'BIOSNAP')

#ht <- get_data('DrugBank')

htnr <- get_data_rmsd('NR')
htic <- get_data_rmsd('IC')
htec <- get_data_rmsd('E')
htgpcr <-  get_data_rmsd('GPCR')
htdrugbank <- get_data_rmsd('DrugBank')
htdavis <- get_data_rmsd('Davis_et_al')
htbindingdb <- get_data_rmsd('BindingDB')
htbiosnap <- get_data_rmsd('BIOSNAP')


require(gridExtra)

file_fig_hist <- 'rmsd_panel.pdf'
pdf(file = file_fig_hist,  width=10, height=6)
  # ggdraw() +
  # draw_plot(htnr, x = 0, y = .5, width = .5, height = .5) +
  # draw_plot(htic, x = .5, y = .5, width = .5, height = .5)
  #ggarrange(htnr, htic, htec, htgpcr,htdrugbank, htdavis, htbindingdb, htbiosnap,
  #           common.legend = FALSE, legend = "bottom")
  grid.arrange(htnr, htic, htec, htgpcr, htdrugbank, htdavis, htbindingdb, htbiosnap,
               nrow=2, ncol=4)
dev.off()



## panel for drug tanimotos


get_data_tanimotos <- function(database_name){
  print(database_name)
  folder_path <- '../Results/data_per_dataset'
  folder_db_path <-  file.path(folder_path, database_name)
  #file_prot_sim <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))
  file_drug_sim <- file.path(folder_path, database_name, paste0('drugs_sim_', database_name, '.csv'))

  # load
  print(file_drug_sim)
  mat_sim_drugs <- read.table(file_drug_sim, sep = ";", header=T)
  mat_sim_drugs <- as.matrix(mat_sim_drugs)

  print( dim(mat_sim_drugs))
  hits_drug = mat_sim_drugs[upper.tri(mat_sim_drugs, diag = FALSE)]
  #
  tit <- database_name
  data_hist <- data.frame(hits_drug)
  # ggplot
  ggplot(data_hist, aes(hits_drug)) +
    ggtitle(tit)+
    labs(x =  "Drug Similarity (Tanimoto score)") +             
    geom_histogram(bins=60, colour="#0c0f0a", fill="#4472C4")+
    #geom_density()+
    theme_classic()+
    theme(plot.title = element_text(hjust = 0.5, size=8), 
          axis.title = element_text(size=5),
          axis.text = element_text(size=5))   
}


#  c('NR', 'IC', 'GPCR', 'E', 'DrugBank', 'Davis_et_al', 'BindingDB', 'BIOSNAP')

dnr <- get_data_tanimotos('NR')
dic <-  get_data_tanimotos('IC')
dgpcr <-  get_data_tanimotos('GPCR')
de <-  get_data_tanimotos('E')
ddrugbank <- get_data_tanimotos('DrugBank')
ddavis <- get_data_tanimotos('Davis_et_al')
dbinding <- get_data_tanimotos('BindingDB')
dbiosnap <- get_data_tanimotos('BIOSNAP')


require(gridExtra)

file_fig_drugs <- 'tani_panel.pdf'
pdf(file = file_fig_drugs,  width=10, height=6)
  grid.arrange(dnr, dic, dgpcr, de, ddrugbank, ddavis, dbinding, dbiosnap,
               nrow=2, ncol=4)
dev.off()

