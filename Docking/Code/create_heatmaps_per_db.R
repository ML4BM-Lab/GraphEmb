
# install.packages("devtools")
#install_github("jokergoo/ComplexHeatmap")
#source("https://bioconductor.org/biocLite.R")
if(!require(ComplexHeatmap)) biocLite("ComplexHeatmap")
if(!require(circlize)) install.packages('circlize')

library(devtools)
library(circlize)
# better rasterization
library(magick)

##################
# Load files
database_name <- 'DrugBank'
print(database_name)
folder_path <- '../Results/data_per_dataset'
folder_db_path <-  file.path(folder_path, database_name)
# files
file_drug_sim <- file.path(folder_path, database_name, paste0('drugs_sim_', database_name, '.csv'))
file_drug_annot <- file.path(folder_path, database_name, paste0('drugs_annot_', database_name, '.csv'))

file_prot_rmsd <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))
file_prot_annot <- file.path(folder_path, database_name, paste0('prots_annot_', database_name, '.csv'))

## LOAD STUFF
# DRUGS
# load drug annots
annots <- read.table(file_drug_annot, header = TRUE, sep = ";", quote = "")
annot_king <- annots$kingdom
annot_superclass <- annots$superclass
annot_class <- annots$class
annot_subclass <- annots$subclass
# load mat sim
mat_sim_drugs <- read.table(file_drug_sim, sep = ";", header=T)
# transform to matrix
mat_sim_drugs <- as.matrix(mat_sim_drugs)

# PROTEINS
# Mat rmsd protein
mat_rmsd_proteins <- read.table(file_prot_rmsd, sep = ";", header=T)
mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
# Protein annotation
annots_prot <- read.table(file_prot_annot, header = TRUE, sep = ";", quote = "")


#### COLOR SETTINGS
library(RColorBrewer)

# color gradient of matix
col_fun <- colorRamp2(c(0, 1), c("#126ac8", "white"))
col_fun(seq(-150, 150))

# For drug classes
# Color for Kingdom
tags_king = unique(annots$kingdom)[unique(annots$kingdom) != "-"]

if (length(tags_king)==1) { 
    color_vec_kig = c('#b80050')
    } else if (length(tags_king)==2) {
        color_vec_kig = c('#b80050','#23074b')
} 
color_kingdom = setNames(color_vec_kig, tags_king)
color_kingdom <- c(color_kingdom, setNames('#ffffff', "-"))

# Load Paletes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# superclass
tags_superclass = unique(annots$superclass)[unique(annots$superclass) != "-"]
color_superclass = setNames(col_vec_super[1:length(tags_superclass)], tags_superclass)
color_superclass <- c(color_superclass, setNames('#FFFFFF', "-"))

# class
tags_class = unique(annots$class)[unique(annots$class) != "-"]
color_class = setNames(col_vec_super[1:length(tags_class)], tags_class)
color_class <- c(color_class, setNames('#FFFFFF', "-"))

# subclass
tags_subclass = unique(annots$subclass)[unique(annots$subclass) != "-"]
color_subclass = setNames(col_vec_super[1:length(tags_subclass)], tags_subclass)
color_subclass <- c(color_subclass, setNames('#FFFFFF', "-"))



# if only kindgom and superclass
# anno_df = data.frame(kingdom = annot_king,
#                     superclass = annot_superclass)

# ha <- HeatmapAnnotation(df = anno_df,
#         simple_anno_size = unit(15, "mm"), width = NULL
#         )

# ------------- All

anno_df = data.frame(kingdom = annot_king,
                    superclass = annot_superclass,
                    class = annot_class,
                    subclass = annot_subclass)

ha <- HeatmapAnnotation(df = anno_df,
        col = list(kingdom = color_kingdom,
                superclass = color_superclass,
                class = color_class,
                subclass = color_subclass),
        simple_anno_size = unit(15, "mm"), width = NULL
        )



file_fig <- paste0(folder_db_path, '/heatmap_drugs_', database_name, '.pdf')
pdf(file = file_fig, width = 20, height = 18)
Heatmap(mat_sim_drugs, col = col_fun, name = 'Drug Similarity',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        top_annotation = ha)
dev.off()




### PROTEIN Matrix

# color gradient of matix
col_fun <- colorRamp2(c(max(mat_rmsd_proteins), 5,0), c("#14487f", "#2586ee","white"))
col_fun(seq(-100, 100))

# Source
tags_source = unique(annots_prot$source)
if (length(tags_source)==1) { 
    color_vec_source = c('#b80050')
    } else if (length(tags_source)==2) {
        color_vec_source = c('#b80050','#23074b')
        }else if (length(tags_source)>2) {
        color_vec_source = col_vec_super[1:length(tags_source)]
}
color_source = setNames(color_vec_source, tags_source)

### Protein class
tag_prot_class = unique(annots_prot$Class)
if (length(tag_prot_class)==1) { 
    color_vec_protclass = c('#b3e065')
    } else if (length(tag_prot_class)==2) {
        color_vec_protclass = c('#b3e065','#9a75cd')
        }else if (length(tag_prot_class)>2) {
        color_vec_protclass = col_vec_super[1:length(tag_prot_class)]
}
color_prot_class= setNames(color_vec_protclass, tag_prot_class)


### only for E 
# brendas = annots_prot$ec
# list_first_b = c()
# for (i in brendas){
#         list_first_b <- append(list_first_b ,substr(as.character(i), 1,1))
# }

# tag_prot_class = unique(list_first_)
# if (length(tag_prot_class)==1) { 
#     color_vec_protclass = c('#b3e065')
#     } else if (length(tag_prot_class)==2) {
#         color_vec_protclass = c('#b3e065','#9a75cd')
#         }else if (length(tag_prot_class)>2) {
#         color_vec_protclass = brewer.pal(n = length(tag_prot_class), name = "Set2")
# }
# color_prot_class= setNames(color_vec_protclass, tag_prot_class)


df_anot_top = data.frame(source = annots_prot$source)

top_annot <- HeatmapAnnotation(df = df_anot_top,
        col = list(source = color_source),
        simple_anno_size = unit(10, "mm"), width = NULL)

#
df_left_top = data.frame(mol_funct = annots_prot$Class)

left_annot <- rowAnnotation(df = df_left_top, lengths = anno_barplot(annots_prot$length, axis = FALSE, 
                                        which = "column", facing='outside'),
        col = list(mol_funct = color_prot_class),
        simple_anno_size = unit(10, "mm"), width = NULL)



file_fig <- paste0(folder_db_path, '/heatmap_prots_', database_name, '.pdf')
pdf(file = file_fig)
Heatmap(mat_rmsd_proteins, col = col_fun, name = 'Protein RMSD',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        top_annotation = top_annot, left_annotation=left_annot
        )
dev.off()




########

print('1- xxxxx 0 xxxx')

# ##############
# # Distribution of Pairwise Similitude

library('ggplot2')
hits = mat_sim_drugs[upper.tri(mat_sim_drugs, diag = FALSE)]

data <- data.frame(hits)
#head(data)

file_fig <- paste0(folder_db_path, '/hist_drugs_', database_name, '.pdf')
pdf(file = file_fig)
ggplot(data, aes(hits)) +
  theme(plot.title = element_text(hjust = 0.5)) +     
  ggtitle("Pairwise similarities histogram")+
  labs(x= "Drug Similarty") +             
  geom_histogram(aes(y = ..density..), bins=80, colour="black", fill="#f68080")+
  geom_density()
dev.off()

# again but histogram rmsd 

hits_prot = mat_rmsd_proteins[upper.tri(mat_rmsd_proteins, diag = FALSE)]
data <- data.frame(hits_prot)
#head(data)

file_fig <- paste0(folder_db_path, '/hist_prots_', database_name, '.pdf')
pdf(file = file_fig)
ggplot(data, aes(hits_prot)) +
  theme(plot.title = element_text(hjust = 0.5)) +     
  ggtitle("Histogram of protein RMSD")+
  labs(x = "RMSD (Å)") +             
  geom_histogram(aes(y = ..density..), bins=60, colour="black", fill="#f9b16e")+
  geom_density()
dev.off()

print('2- xxxxx 0 xxxx')
