
if(!require(ComplexHeatmap)) biocLite("ComplexHeatmap")
#if(!require(circlize)) install.packages('circlize')

library(devtools)
library(circlize)
# better rasterization
library(magick)
library(ggplot2)
library(RColorBrewer)


#############
# DRUGS
# Paths
folder_path <- '../Results/data_per_dataset/full_data'
FILE_ANNOT <- file.path(folder_path, 'df_annot_drugs.csv')
FILE_SIM <- file.path(folder_path, 'tani_drugs.csv')

# Load  
annots_drugs <- read.table(FILE_ANNOT, header = TRUE, sep = ";", quote = "")
mat_drug_sim <- read.table(FILE_SIM, sep = ";", header=T)
mat_drug_sim <- as.matrix(mat_drug_sim)


########################
### PLOT HISTOGRAM
# First generate histogram
hits_drug = mat_drug_sim[upper.tri(mat_drug_sim, diag = FALSE)]
data <- data.frame(hits_drug)
#head(data)

file_fig_hist <- file.path(folder_path, 'hist_drugs.pdf')

pdf(file = file_fig_hist)
ggplot(data, aes(hits_drug)) +
  ggtitle("Histogram of drug pairwise similarity")+
  labs(x = "Drug Similarity (Tanimoto score)") +             
  geom_histogram(aes(y = ..density..), bins=60, colour="#0c0f0a", fill="#d04088")+
  geom_density()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))   
dev.off()


### Generate Heatmap

############### COLORS
# color gradient of matix
#pale <- c("#03045e","#0077b6","#00b4d8","#90e0ef","#caf0f8"])
#col_fun <- colorRamp2(c(0, 0.25, 0.5,0.75,1), pale)

pale <- c("#133c55","#386fa4","#78b7e2","#f5f5f5")
col_fun <- colorRamp2(c(0, 0.5,0.75,1), pale)
col_fun(seq(-150, 150))

# col_fun = circlize::colorRamp2(seq(0, 1, by=(0-(-1))/(150-1)),
#  , colors=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 11, name = "Blues")))(150))


### Annotations
annot_kingdom <- annots_drugs$kingdom
tags_king <- unique(annot_kingdom)[unique(annot_kingdom) != ""] # change again for -

if (length(tags_king)==1) { 
    color_vec_kig = c('#17b46a')
    } else if (length(tags_king)==2) {
        color_vec_kig = c('#17b46a','#ee451b')
} 
color_kingdom <- setNames(color_vec_kig, tags_king)
color_kingdom <- c(color_kingdom, setNames('#ffffff', ""))

# Load Paletes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# superclass
tags_superclass = unique(annot_superclass)[unique(annot_superclass) != "-"]
color_superclass = setNames(col_vec_super[1:length(tags_superclass)], tags_superclass)
color_superclass <- c(color_superclass, setNames('#FFFFFF', "-"))


######### Molecular Descriptors

names(annots_drugs)

logp <- annots_drugs$log_p
numhet <- annots_drugs$numhet
molwt <- annots_drugs$molwt

pale_molwt <- c("#78f47c", "#f5f5f5")

col_fun_molwt <- colorRamp2(c(max(molwt), mean(molwt), min(molwt)), c('black', 'orange','white'))
col_fun_molwt(seq(-200, 200))


col_fun_logp <- colorRamp2(c(max(logp), mean(logp), min(logp)), c('black', '#589cdb','white'))
col_fun_logp(seq(-200, 200))


col_fun_numhet <- colorRamp2(c(max(numhet), mean(numhet), min(numhet)), c('black', '#b740ad','white'))
col_fun_numhet(seq(-200, 200))

# create dataframe for annotations
#if only kindgom and superclass
anno_df = data.frame(kingdom = annot_kingdom[1:400],
                     superclass = annot_superclass[1:400])


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

anno_df_moldes = data.frame(molwt = molwt[1:400],
                            logp = logp[1:400],
                            numhet = numhet[1:400])

left_ha <- rowAnnotation(df=anno_df_moldes, 
                        col= list(molwt= col_fun_molwt,
                                  logp = col_fun_logp,
                                   numhet = col_fun_numhet), 
                        simple_anno_size = unit(15, "mm"))


pdf(file = "../Results/test_hmap_v3.pdf", width = 20, height = 18)
Heatmap(mat_drug_sim[1:400, 1:400], col = col_fun, name = 'Drug Similarity',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE, 
        heatmap_legend_param = list(at = c(0, 0.25, 0.5,0.75,1),
                                legend_direction = "vertical",
                                legend_height = unit(20, "cm") 
                                ),
        top_annotation = ha, 
        right_annotation = row_ha,
        left_annotation = left_ha,
        use_raster = TRUE
)
dev.off()


