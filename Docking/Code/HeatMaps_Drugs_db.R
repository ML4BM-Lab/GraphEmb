
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
database_name <- 'BindingDB'

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
# hits_drug = mat_sim_drugs[upper.tri(mat_sim_drugs, diag = FALSE)]
# data <- data.frame(hits_drug)

# file_fig <- paste0(folder_db_path, '/hist_drugs_', database_name, '.pdf')

# pdf(file = file_fig)
# ggplot(data, aes(hits_drug)) +
#   ggtitle(paste0("Dataset: ", database_name))+
#   labs(x = "Drug Similarity (Tanimoto score)") +             
#   geom_histogram(aes(y = ..density..), bins=60, colour="#0c0f0a", fill="#d04088")+
#   geom_density()+
#   theme_classic()+
#   theme(plot.title = element_text(hjust = 0.5))   
# dev.off()


### HEATMAP 

# color heatmap
pale <- c("#133c55","#386fa4","#78b7e2","#f5f5f5")
col_fun <- colorRamp2(c(0, 0.5,0.75,1), pale)
col_fun(seq(-300, 300))



### Annotations

# Not using annot_kingdom
# tags_king <- unique(annot_kingdom)[unique(annot_kingdom) != "-"]

# if (length(tags_king)==1) { 
#     color_vec_kig = c('#17b46a')
#     } else if (length(tags_king)==2) {
#         color_vec_kig = c('#17b46a','#ee451b')
# } 
# color_kingdom <- setNames(color_vec_kig, tags_king)
# color_kingdom <- c(color_kingdom, setNames('#ffffff', "-"))


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
## Not using subclass
# tags_subclass = unique(annot_subclass)[unique(annot_subclass) != "-"]
# color_subclass = setNames(col_vec_super[1:length(tags_subclass)], tags_subclass)
# color_subclass <- c(color_subclass, setNames('#FFFFFF', "-"))


#### Molecular Descriptors


logp <- annots_drugs$log_p
numhet <- annots_drugs$numhet
molwt <- annots_drugs$molwt

pale_molwt <- c("#78f47c", "#f5f5f5")

col_fun_molwt <- colorRamp2(c(max(molwt), mean(molwt)*2, mean(molwt), mean(molwt)/2, min(molwt)), 
                                c("#cc5803","#e2711d","#ff9505","#ffb627","#ffc971"))
col_fun_molwt(seq(-200, 200))

col_fun_logp <- colorRamp2(c(max(logp), mean(logp), min(logp)), c('black', '#589cdb','white'))
col_fun_logp(seq(-200, 200))
col_fun_numhet <- colorRamp2(c(max(numhet), mean(numhet), min(numhet)), c('black', '#b740ad','white'))
col_fun_numhet(seq(-200, 200))



# create dataframe for annotations
#if only kindgom and superclass

anno_df = data.frame(kingdom = annot_kingdom,
                     superclass = annot_superclass,
                     class = annot_class,
                     subclass = annot_subclass)


ha <- HeatmapAnnotation(df = anno_df,
                        col = list(kingdom = color_kingdom,
                                superclass = color_superclass,
                                class = color_class,
                                subclass = color_subclass),
                        simple_anno_size = unit(15, "mm"), 
                        width = NULL
                        )


row_ha = rowAnnotation(df = anno_df, 
                        col = list(kingdom = color_kingdom,
                                superclass = color_superclass,
                                class = color_class,
                                subclass = color_subclass),
                        simple_anno_size = unit(15, "mm"),
                        show_legend = FALSE,
                        show_annotation_name = FALSE)


anno_df_moldes = data.frame(molwt = molwt,
                            logp = logp,
                            numhet = numhet)

left_ha <- rowAnnotation(df=anno_df_moldes, 
                        col= list(molwt= col_fun_molwt,
                                  logp = col_fun_logp,
                                   numhet = col_fun_numhet), 
                        simple_anno_size = unit(15, "mm"))




file_fig <- paste0(folder_db_path, '/heatmap_drugs_', database_name, '_new.png')
png(file = file_fig, width = 2000, height = 2000)
Heatmap(mat_sim_drugs, col = col_fun, name = 'Drug Similarity',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE, 
        heatmap_legend_param = list(at = c(0, 0.25, 0.5,0.75,1),
                                legend_direction = "vertical",
                                legend_height = unit(5, "cm") 
                                ),
        top_annotation = ha, 
        left_annotation = row_ha,
        right_annotation = left_ha,
        use_raster = FALSE
)
dev.off()


print('xxxxx')
print(database_name)
