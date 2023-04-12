##
if(!require(ComplexHeatmap)) biocLite("ComplexHeatmap")
#if(!require(circlize)) install.packages('circlize')

library(devtools)
library(circlize)
# better rasterization
library(magick)
library(ggplot2)
library(RColorBrewer)
library(plyr)


##################
# LOAD FILES
database_name <- 'DrugBank' # 'E' 'IC', 'GPCR','NR', 'Davis_et_al'

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


###  

# Color heatmap
pale <- c("#133c55","#386fa4","#78b7e2","#f5f5f5")
col_fun <- colorRamp2(c(0, 0.5,0.75,1), pale)
col_fun(seq(-300, 300))

# Load Paletes
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

godsnot_64 = c(
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
    "#997D87", "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF",
    "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
    "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
    "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
    "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
    "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
    "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
    "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
    "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
    "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
    "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
    "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
    "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

col_vec_super <- godsnot_64

thres_drugs <- 30 # for <= E
counter_superclass <- count(annot_superclass)
selected_superclass <- counter_superclass["freq"] >= thres_drugs # from plyr

which(selected_superclass)

selected_annots <- counter_superclass[which(selected_superclass),'x']

a_sup <- annot_superclass # copy

#print(is.element(annot, selected_annots)

a_sup[!a_sup %in% selected_annots] <- 'Other'

a_sup[which(a_sup=="-")] <- "NotAnnotated"

#as.data.frame(table(annot_superclass))

# superclass
tags_superclass = unique(annot_superclass)[(unique(annot_superclass) != "NotAnnotated") & (unique(annot_superclass) != "Other")]

color_superclass <- setNames(col_vec_super[1:length(tags_superclass)], tags_superclass)
color_superclass <- c(color_superclass, setNames('#808080', "NotAnnotated"), setNames('#000000', "Other"))


# class
# Missing class other here
# Not class for E
# count(tags_class)
# counter_class <- count(tags_class)
# selected_class <- counter_class["freq"] > thres_drugs # from plyr


# tags_class = unique(annot_class)[unique(annot_class) != "-"]
# color_class = setNames(col_vec_super[1+6:(length(tags_class)+5)], tags_class)
# color_class <- c(color_class, setNames('#FFFFFF', "-"))



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

anno_df = data.frame(
                     superclass = a_sup
                     #,class = annot_class
                     )

anno_df_moldes = data.frame(molwt = molwt,
                            logp = logp,
                            numhet = numhet)



## annotations

ha <- HeatmapAnnotation(df = anno_df,
                        col = list(
                                superclass = color_superclass
                                #,class = color_class
                                ),
                        simple_anno_size = unit(15, "mm"), 
                        annotation_name_gp= gpar(fontsize = 25), # name annotation in rows outiside plot
                        width = NULL,
                        annotation_legend_param = list(title = "Superclass",title_position = "topcenter",
                                                        title_gp = gpar(fontsize = 24),
                                                         grid_height = unit(15, "mm"), grid_width=  unit(15, "mm"),
                                                         labels_gp = gpar(col = "black", fontsize = 18)  ) )



row_ha = rowAnnotation(df = anno_df, 
                        col = list(
                                superclass = color_superclass
                                #,class = color_class
                                ),
                        simple_anno_size = unit(15, "mm"),
                        show_legend = FALSE,
                        show_annotation_name = FALSE)



left_ha <- rowAnnotation(df=anno_df_moldes, 
                        col= list(molwt= col_fun_molwt,
                                  logp = col_fun_logp,
                                   numhet = col_fun_numhet), 
                                 annotation_name_gp= gpar(fontsize = 25), # name annotation in rows outiside plot
                                simple_anno_size = unit(15, "mm"),
                                annotation_legend_param = list(legend_height = unit(10, "cm"), legend_direction = "vertical",
                                                                title_position = "lefttop-rot", labels_gp = gpar(fontsize = 18),
                                                                title_gp = gpar(fontsize = 25))
                                )





hmp <- Heatmap(mat_sim_drugs, col = col_fun, name = 'Drug Similarity',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE, 
        heatmap_legend_param = list(at = c(0, 0.25, 0.5,0.75,1),
                                title_position = "lefttop-rot", # girar titutlo leyenda
                                title_gp = gpar(fontsize = 25), # tamaño titulo leyenda
                                legend_direction = "vertical",
                                legend_height = unit(10, "cm"),
                                labels_gp = gpar(fontsize = 18) # numbers in legend bar, but not title
                                ),
        #row_names_gp = gpar(fontsize = 20),
        #row_names_max_width = unit(6, "cm"),
        top_annotation = ha, 
        left_annotation = row_ha,
        right_annotation = left_ha,
        use_raster = FALSE
        )



file_lfig <- paste0(folder_db_path, '/heatmap_drugs_', database_name, '_new.png')

png(file = file_fig, 
                width     = 20,
                height    = 12,
                units     = "in",
                res       = 400,
                pointsize = 4)
hmp
dev.off()

#width = 2400, height = 2000, res= 330) # 2000


#file_fig_pdf = paste0(folder_db_path, '/heatmap_drugs_', database_name, '_new.pdf')
#pdf(file, width=10, height=10)

print('xxxxx')
print(database_name)
