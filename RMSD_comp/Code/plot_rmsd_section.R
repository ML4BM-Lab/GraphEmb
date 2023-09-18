library(ggplot2)



database_name <- 'BIOSNAP'
print(database_name)

folder_path <- '../Results/data_per_dataset'
folder_db_path <-  file.path(folder_path, database_name)
file_prot_sim <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))

# load
print(file_prot_sim)
mat_rmsd_proteins <- read.table(file_prot_sim, sep = ";", header=T)
mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
print( dim(mat_rmsd_proteins))

hits_prot = mat_rmsd_proteins[upper.tri(mat_rmsd_proteins, diag = FALSE)]
data_hist <- data.frame(hits_prot)

file_fig_hist <- paste0(folder_db_path, '/hist_prot_', database_name, '.pdf')

file_fig_hist <- paste0('test_plot', database_name, '.pdf')


pdf(file = file_fig_hist, width=12, height=8)
ggplot(data_hist, aes(hits_prot)) +
  ggtitle("Histogram of RMSD distribution")+
  labs(x = "RMSD (Å)") +             
  geom_histogram(bins=160, colour="#0c0f0a", fill="#f18f01", alpha=0.95, size=0.01)+
#  geom_density()+
  #xlim(0,40)+
  scale_x_continuous(limits = c(0,22),
                     breaks = seq(0,22,2))+
  #scale_y_continuous(breaks=seq(0,100000,5000))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size=20))   
dev.off()


data_mol_biosnap <- read.csv(file = '../Results/moltrans_biosnap.csv', sep=';')
df <-data.frame(data_mol_biosnap)

df[,3]
plot(x, y1, type = "b", pch = 19, 
     col = "red", xlab = "x", ylab = "y")