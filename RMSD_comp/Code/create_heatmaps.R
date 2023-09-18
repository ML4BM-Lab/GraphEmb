
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
FILE_ANNOT <- "../Data/df_annot_drugbank.csv"
FILE_SIM <- "../Data/tani_drugbank.csv"
# tests
#FILE_SIM<- "../Data/data_per_dataset/NR/drugs_NR.csv"
# Load annotations (drug classification)
annots <- read.table(FILE_ANNOT, header = TRUE, sep = ";", quote = "")
annot_king <- annots$kingdom
annot_superclass <- annots$superclass
annot_class <- annots$class
annot_subclass <- annots$subclass


mat_sim <- read.table(FILE_SIM, sep = ";", header=T)


# colnames(mat_sim) <- sub("^X", "", colnames(mat_sim))
# rownames(mat_sim) <- colnames(mat_sim)

## sizes for test! 
# ext_test = 300
# mat_sim <- mat_sim[1:ext_test,1:ext_test]
# annot_king <- annots$kingdom[1:ext_test]
# annot_superclass <- annots$superclass[1:ext_test]
# annot_class <- annots$class[1:ext_test]
# annot_subclass <- annots$subclass[1:ext_test] # a lot 


mat_sim <- as.matrix(mat_sim)



#### COLOR SETTINGS

# color gradient of matix
col_fun <- colorRamp2(c(0, 1), c("#187fed", "white"))
col_fun(seq(-150, 150))

# Color for Kingdom
tags_king = unique(annots$kingdom)[unique(annots$kingdom) != "-"]
color_vec_kig = c('#b80050', '#23074b')
color_kingdom = setNames(color_vec_kig, tags_king)
color_kingdom <- c(color_kingdom, setNames('#ffffff', "-"))


library(RColorBrewer)
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
# superclass
tags_superclass = unique(annots$superclass)[unique(annots$superclass) != "-"]
color_superclass = setNames(col_vec_super[1:length(tags_superclass)], tags_superclass)
color_superclass <- c(color_superclass, setNames('#FFFFFF', "-"))

# subclass
tags_class = unique(annots$class)[unique(annots$class) != "-"]
color_class = setNames(col_vec_super[1:length(tags_class)], tags_class)
color_class <- c(color_class, setNames('#FFFFFF', "-"))


# if only kindgom and superclass
anno_df = data.frame(kingdom = annot_king,
                    superclass = annot_superclass)

ha <- HeatmapAnnotation(df = anno_df,
        col = list(kingdom = color_kingdom,
                  superclass = color_superclass),
        simple_anno_size = unit(15, "mm"), width = NULL
        )


# ------------- All
# anno_df = data.frame(kingdom = annot_king,
#                     superclass = annot_superclass,
#                     class = annot_class)

# ha <- HeatmapAnnotation(df = anno_df,
#         col = list(kingdom = color_kingdom,
#                   superclass = color_superclass,
#                   class = color_class),
#         simple_anno_size = unit(15, "mm"), width = NULL
#         )



pdf(file = "Rplots_full_V3d.pdf", width = 20, height = 18)
Heatmap(mat_sim, col = col_fun, name = 'Drug Similarity',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        top_annotation = ha)

dev.off()


print('xxxxx 0 xxxx')

# ##############
# # Distribution of Pairwise Similitude

# library('ggplot2')
# hits = mat_sim[upper.tri(mat_sim, diag = FALSE)]

# data <- data.frame(hits)
# #head(data)
# pdf(file = "full_hist_drugs.pdf")
# ggplot(data, aes(hits)) +
#   theme(plot.title = element_text(hjust = 0.5)) +     
#   ggtitle("Pairwise similarities histogram")+
#   labs(x= "Drug Similarty") +             
#   geom_histogram(aes(y = ..density..), bins=80, colour="black", fill="#d4a6ed")+
#   geom_density()
# dev.off()


file_prot_rmsd <- '../Data/full_rmsd.csv'
file_prot_annot <- '../Data/full_annotated_proteins.csv'



mat_sim_rmsd <- read.table(file_prot_rmsd, sep = ";", header=T)
annots_prot <- read.table(file_prot_annot, header = TRUE, sep = ";", quote = "")

mat_sim_rmsd <- as.matrix((mat_sim_rmsd))

#### COLOR SETTINGS
library(RColorBrewer)

# color gradient of matix
col_fun <- colorRamp2(c(max(mat_sim_rmsd), 50, 10, 5, 0), c("#07090b", "#122e4a","#649fde","#6ee7dd","white"))
col_fun(seq(-200, 200))

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

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



file_fig <- 'test_all_prot_v2.pdf'
pdf(file = file_fig)
Heatmap(mat_sim_rmsd, col = col_fun, name = 'Protein RMSD',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        top_annotation = top_annot, left_annotation=left_annot
        )
dev.off()


library('ggplot2')

hits_prot = mat_sim_rmsd[upper.tri(mat_sim_rmsd, diag = FALSE)]
data <- data.frame(hits_prot)
#head(data)

file_fig <- 'test_hist_prot.pdf'
pdf(file = file_fig)
ggplot(data, aes(hits_prot)) +
  theme(plot.title = element_text(hjust = 0.5)) +     
  ggtitle("Histogram of protein RMSD")+
  labs(x = "RMSD (Å)") +             
  geom_histogram(aes(y = ..density..), bins=60, colour="black", fill="#f9b16e")+
  geom_density()
dev.off()