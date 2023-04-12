
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
database_name <- 'DrugBank' #E, GPCR, IC, NR, Davis_et_al,BindingDB

print(database_name)
folder_path <- '../Results/data_per_dataset'
folder_db_path <-  file.path(folder_path, database_name)

# FILES
file_prot_sim <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))
file_prot_annot <- file.path(folder_path, database_name, paste0('prots_annot_', database_name, '.csv'))


# load
mat_rmsd_proteins <- read.table(file_prot_sim, sep = ";", header=T)
mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
dim(mat_rmsd_proteins)
# Protein annotation
annots_prot <- read.table(file_prot_annot, header = TRUE, sep = ";", quote = "")
names(annots_prot)




####### HEATMAP

# color gradient of matix
if (max(mat_rmsd_proteins)>50){
    col_fun <- colorRamp2(c(max(mat_rmsd_proteins), 50, 10, 5, 0), c("#07090b", "#122e4a","#649fde","#6ee7dd","white"))
    col_fun(seq(-200, 200))
    list_legend = c(max(mat_rmsd_proteins), 50, 10, 5, 0)}else{
    col_fun <- colorRamp2(c(max(mat_rmsd_proteins), 10, 0), c("#07090b", "#438ddc","white"))
    col_fun(seq(-200, 200))
    list_legend = c(max(mat_rmsd_proteins), 10, 0)
    }


# Load Paletes
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))



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




# SOURCE 
tags_source = unique(annots_prot$source)
if (length(tags_source)==1) { 
    color_vec_source = c('#b80050')
    } else if (length(tags_source)==2) {
        color_vec_source = c('#b80050','#23074b')
        }else if (length(tags_source)>2) {
        color_vec_source = col_vec_super[1:length(tags_source)]
}

color_source = setNames(color_vec_source, tags_source)


### Protein Molecular Function

# for enzimes not

thres_mf <- 10 # for largets datasets
counter_molfunc <- count(annots_prot$mol_func)

selected__molfunc <- counter_molfunc["freq"] >= thres_mf # from plyr
which(selected__molfunc)

selected_molf <- counter_molfunc[which(selected__molfunc),'x']

a_mf <- annots_prot$mol_fun # copy

a_mf[!a_mf %in% selected_molf] <- 'other'


tags_mf = unique(a_mf)[(unique(a_mf) != "notAnnotated") & (unique(a_mf) != "other")]

color_prot_class <- setNames(col_vec_super[1:length(tags_mf)], tags_mf)
color_prot_class <- c(color_prot_class, setNames('#808080', "notAnnotated"), setNames('#000000', "other"))


# tag_prot_class = unique(anno_df_moldes)

# if (length(tag_prot_class)==1) { 
#     color_vec_protclass = c('#b3e065')
#     } else if (length(tag_prot_class)==2) {
#         color_vec_protclass = c('#b3e065','#9a75cd')
#         }else if (length(tag_prot_class)>2) {
#         color_vec_protclass = col_vec_super[1+9:(length(tag_prot_class)+8)]
# }

# color_prot_class= setNames(color_vec_protclass, tag_prot_class)

#color_prot_class <- c(color_prot_class, setNames('#ffffff', "-"))



# Enzyme code
ec_code <- annots_prot$EC
tags_ec <- sort(unique(ec_code)[unique(ec_code) != "-"]) # change again for -
if (length(tags_ec)!=0){
    color_list_ec = c('#ea481f', '#ffad0a', '#6dcd3c', '#f50093', '#009eb4', '#6a1145', '#0f27bf')
    color_ec <- setNames(color_list_ec[1:length(tags_ec)], tags_ec)
    color_ec <- c(color_ec, setNames('#000000', "-"))
}


#### Create annotations



## TOP ANNOTATION
# SOURCE AND MOLECULAR FUNCTION
annot_top = data.frame(Source = annots_prot$source,
                       Function  = a_mf) # changed here !! instead of annots_prot$annot_ annots_prot$mol_func

ha_top <- HeatmapAnnotation(df = annot_top,
        col = list(Source = color_source,
                    Function = color_prot_class),
        show_legend = c(TRUE, FALSE),
        width = NULL,
        annotation_name_gp= gpar(fontsize = 25), # name annotation in rows outiside plot
        simple_anno_size = unit(15, "mm"),
        annotation_legend_param = list(legend_height = unit(10, "cm"), legend_direction = "vertical",
                                        title_position = "lefttop-rot", labels_gp = gpar(fontsize = 18),
                                        title_gp = gpar(fontsize = 20))
                                        )


# LENGTH AND FUNCTION
# min(annots_prot$length) changed for 0
# max(annots_prot$length) changed for 1000
lengths <- as.double(annots_prot$length)

col_fun_len <- colorRamp2(c(max(annots_prot$length), 800, 400,200, min(annots_prot$length)), 
c('#3d1e2e','#5d2e46', '#bc928a','#d0ada7', '#eedfd6'))
col_fun_len(seq(-150, 150))




#### ANNOTATION LEFT
### for enzyme code 
annot_left = data.frame(Length = lengths,
                        Function  = a_mf,
                        EC = annots_prot$EC
                        )


left_annot <- rowAnnotation(df = annot_left,
                col = list(Length = col_fun_len,
                         Function = color_prot_class,
                         EC = color_ec),
                         width = NULL,
                         show_annotation_name=TRUE,
                          annotation_name_gp= gpar(fontsize = 25), # name annotation in rows outiside plot
                          simple_anno_size = unit(15, "mm"),
                          annotation_legend_param = list(legend_height = unit(5, "cm"), legend_direction = "vertical",
                                                          title_position = "lefttop-rot", labels_gp = gpar(fontsize = 18),
                                                          title_gp = gpar(fontsize = 20))
                                                          )


#No enzyme code                
# annot_left = data.frame(Length = lengths,
#                         Function  = a_mf
#                         )

# left_annot <- rowAnnotation(df = annot_left,
#                 col = list(Length = col_fun_len,
#                         Function = color_prot_class),
#                 width = NULL,
#                 show_annotation_name=TRUE,
#                 annotation_name_gp= gpar(fontsize = 25), # name annotation in rows outiside plot
#                 simple_anno_size = unit(15, "mm"),
#                 annotation_legend_param = list(legend_height = unit(5, "cm"), legend_direction = "vertical",
#                                                 title_position = "lefttop-rot", labels_gp = gpar(fontsize = 18),
#                                                 title_gp = gpar(fontsize = 20))
#                                                 )


### SAVE

#file_fig_heat <- file.path(folder_path, 'heatmap_prots.png')


hmp <- Heatmap(mat_rmsd_proteins, 
        col = col_fun, name = 'Protein RMSD',
        clustering_distance_rows = "euclidean",
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        heatmap_legend_param = list(at = list_legend,
                                title_position = "lefttop-rot", # girar titutlo leyenda
                                title_gp = gpar(fontsize = 20), # tamaño titulo leyenda
                                legend_direction = "vertical",
                                legend_height = unit(8, "cm"),
                                labels_gp = gpar(fontsize = 18) # numbers in legend bar, but not title
                                ),
        border = TRUE,
        top_annot = ha_top,
        left_annotation = left_annot,
        )


file_fig_heat <- paste0(folder_db_path, '/heatmap_prots_', database_name, '_new.png')

png(file = file_fig_heat, 
                width     = 20,
                height    = 15,
                units     = "in",
                res       = 400,
                pointsize = 4)
hmp
dev.off()

print('xxxxx')
print(database_name)
