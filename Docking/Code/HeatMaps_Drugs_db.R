
if(!require(ComplexHeatmap)) biocLite("ComplexHeatmap")
#if(!require(circlize)) install.packages('circlize')

library(devtools)
library(circlize)
# better rasterization
library(magick)
library(ggplot2)
library(RColorBrewer)


##################
# LOAD FILES
database_name <- 'E'

print(database_name)
folder_path <- '../Results/data_per_dataset'
folder_db_path <-  file.path(folder_path, database_name)

# FILES
file_drug_sim <- file.path(folder_path, database_name, paste0('drugs_sim_', database_name, '.csv'))
file_drug_annot <- file.path(folder_path, database_name, paste0('drugs_annot_', database_name, '.csv'))


## load 
# Annotations
annots_drugs <- read.table(file_drug_annot, header = TRUE, sep = ";", quote = "")
annot_kingdom <- annots_drugs$kingdom
annot_superclass <- annots_drugs$superclass
annot_class <- annots_drugs$class
annot_subclass <- annots_drugs$subclass

# Load similitude matrix
mat_sim_drugs <- read.table(file_drug_sim, sep = ";", header=T)
mat_sim_drugs <- as.matrix(mat_sim_drugs)


###  HISTOGRAM
hits_drug = mat_sim_drugs[upper.tri(mat_sim_drugs, diag = FALSE)]
data <- data.frame(hits_drug)

file_fig <- paste0(folder_db_path, '/hist_drugs_', database_name, '_v2.pdf')

pdf(file = file_fig)
ggplot(data, aes(hits_drug)) +
  ggtitle(paste0("Dataset: ", database_name))+
  labs(x = "Drug Similarity (Tanimoto score)") +             
  geom_histogram(aes(y = ..density..), bins=60, colour="#0c0f0a", fill="#d04088")+
  geom_density()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))   
dev.off()


### HEATMAP 

# color heatmap
pale <- c("#133c55","#386fa4","#78b7e2","#f5f5f5")
col_fun <- colorRamp2(c(0, 0.5,0.75,1), pale)
col_fun(seq(-150, 150))

### Annotations
tags_king <- unique(annot_kingdom)[unique(annot_kingdom) != "-"]

if (length(tags_king)==1) { 
    color_vec_kig = c('#17b46a')
    } else if (length(tags_king)==2) {
        color_vec_kig = c('#17b46a','#ee451b')
} 
color_kingdom <- setNames(color_vec_kig, tags_king)
color_kingdom <- c(color_kingdom, setNames('#ffffff', "-"))

# Load Paletes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# superclass
tags_superclass = unique(annot_superclass)[unique(annot_superclass) != "-"]
color_superclass = setNames(col_vec_super[1:length(tags_superclass)], tags_superclass)
color_superclass <- c(color_superclass, setNames('#FFFFFF', "-"))

# class
tags_class = unique(annot_class)[unique(annot_class) != "-"]
color_class = setNames(col_vec_super[1+6:(length(tags_class)+5)], tags_class)
color_class <- c(color_class, setNames('#FFFFFF', "-"))

# subclass
tags_subclass = unique(annot_subclass)[unique(annot_subclass) != "-"]
color_subclass = setNames(col_vec_super[1:length(tags_subclass)], tags_subclass)
color_subclass <- c(color_subclass, setNames('#FFFFFF', "-"))



# create dataframe for annotations
#if only kindgom and superclass
anno_df = data.frame(kingdom = annot_kingdom,
                     superclass = annot_superclass)


ha <- HeatmapAnnotation(df = anno_df,
                        col = list(kingdom = color_kingdom,
                                superclass = color_superclass),
                        simple_anno_size = unit(15, "mm"), 
                        width = NULL
                        )


row_ha = rowAnnotation(df = anno_df, 
                        col = list(kingdom = color_kingdom,
                                superclass = color_superclass),
                        simple_anno_size = unit(15, "mm"),
                        show_legend = FALSE,
                        show_annotation_name = FALSE)



anno_df_2 = data.frame(class = annot_class,
                     subclass = annot_subclass)

row_ha_right = rowAnnotation(df = anno_df_2, 
                        col = list(class = color_class,
                                subclass = color_subclass),
                        simple_anno_size = unit(15, "mm"),
                        show_legend = TRUE,
                        show_annotation_name = TRUE)




file_fig <- paste0(folder_db_path, '/heatmap_drugs_', database_name, '_v2.pdf')
pdf(file = file_fig, width = 20, height = 18)
Heatmap(mat_sim_drugs, col = col_fun, name = 'Drug Similarity',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE, 
        heatmap_legend_param = list(at = c(0, 0.25, 0.5,0.75,1),
                                legend_direction = "vertical",
                                legend_height = unit(20, "cm") 
                                ),
        top_annotation = ha, 
        left_annotation = row_ha,
        right_annotation = row_ha_right,
        use_raster = TRUE
)
dev.off()


