
if(!require(ComplexHeatmap)) biocLite("ComplexHeatmap")
#if(!require(circlize)) install.packages('circlize')

library(devtools)
library(circlize)
# better rasterization
library(magick)
library(ggplot2)
library(RColorBrewer)



# Proteins


get_matrix_proteins <- function(database_name){
    print(database_name)
    folder_path <- '../Results/data_per_dataset'
    folder_db_path <-  file.path(folder_path, database_name)
    file_prot_sim <- file.path(folder_path, database_name, paste0('prots_rmsd_', database_name, '.csv'))
    mat_rmsd_proteins <- read.table(file_prot_sim, sep = ";", header=T)
    mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
    mat_rmsd_proteins
    }

get_annot_proteins <- function(database_name){
    folder_path <- '../Results/data_per_dataset'
    folder_db_path <-  file.path(folder_path, database_name)
    file_prot_annot <- file.path(folder_path, database_name, paste0('prots_annot_', database_name, '.csv'))
    annots_prot <- read.table(file_prot_annot, header = TRUE, sep = ";", quote = "")
}


get_list_legend <- function(mat_rmsd_proteins){
    # color gradient of matix
    if (max(mat_rmsd_proteins)>50){
        col_fun <- colorRamp2(c(max(mat_rmsd_proteins), 50, 10, 5, 0), c("#07090b", "#122e4a","#649fde","#6ee7dd","white"))
        col_fun(seq(-200, 200))
        list_legend = c(max(mat_rmsd_proteins), 50, 10, 5, 0)}else{
        col_fun <- colorRamp2(c(max(mat_rmsd_proteins), 10, 0), c("#07090b", "#438ddc","white"))
        col_fun(seq(-200, 200))
        list_legend = c(max(mat_rmsd_proteins), 10, 0)
        }
        list_legend
}

# SOURCE 
get_color_source <- function(annots_prot){
    tags_source = unique(annots_prot$source)
    if (length(tags_source)==1) { 
        color_vec_source = c('#b80050')
        } else if (length(tags_source)==2) {
            color_vec_source = c('#b80050','#23074b')
            }else if (length(tags_source)>2) {
            color_vec_source = col_vec_super[1:length(tags_source)]
    }
    color_source = setNames(color_vec_source, tags_source)
}

get_prot_class <- function(annots_prot){
    ### Protein Molecular Function
    tag_prot_class = unique(annots_prot$mol_func)
    if (length(tag_prot_class)==1) { 
        color_vec_protclass = c('#b3e065')
        } else if (length(tag_prot_class)==2) {
            color_vec_protclass = c('#b3e065','#9a75cd')
            }else if (length(tag_prot_class)>2) {
            color_vec_protclass = col_vec_super[1+9:(length(tag_prot_class)+8)]
    }
    color_prot_class= setNames(color_vec_protclass, tag_prot_class)
}


#################################################
#################################################

list_datasets = c('NR', 'IC', 'GPCR', 'E', 'DrugBank', 'Davis_et_al', 'BindingDB', 'BIOSNAP')

# for (dataset in list_datasets){
#     database_name = dataset
#     #print(database_name)


database_name = 'BIOSNAP'


mat_rmsd_proteins <- get_matrix_proteins(database_name)
annots_prot <- get_annot_proteins(database_name)

#
# Load Paletes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# color heatmap
list_legend <- get_list_legend(mat_rmsd_proteins)

# colors annotations
color_source <- get_color_source(annots_prot)
color_prot_class <- get_prot_class(annots_prot)


## Enzyme code if there are enzymes
ec_code <- annots_prot$EC
tags_ec <- sort(unique(ec_code)[unique(ec_code) != "-"]) # change again for -
if (length(tags_ec)!=0){
    color_list_ec = c('#ea481f', '#ffad0a', '#6dcd3c', '#f50093', '#009eb4', '#6a1145', '#041f1e')
    color_ec <- setNames(color_list_ec[1:length(tags_ec)], tags_ec)
    color_ec <- c(color_ec, setNames('#ffffff', "-"))
}

# LENGTH AND FUNCTION
# min(annots_prot$length) changed for 0
# max(annots_prot$length) changed for 1000
lengths <- as.double(annots_prot$length)

col_fun_len <- colorRamp2(c(max(annots_prot$length), 800, 400,200, min(annots_prot$length)), 
                            c('#3d1e2e','#5d2e46', '#bc928a','#d0ada7', '#eedfd6'))
col_fun_len(seq(-150, 150))


# if there are encymes, annotation
if (length(tags_ec)!=0){
    annot_left = data.frame(Length = lengths,
                            Function  = annots_prot$mol_func,
                            EC = annots_prot$EC
                            )

    left_annot <- rowAnnotation(df = annot_left,
                    col = list(Length = col_fun_len,
                            Function = color_prot_class,
                            EC = color_ec),
                    simple_anno_size = unit(10, "mm"), 
                    width = NULL,
                    show_annotation_name=TRUE)}else{
        annot_left = data.frame(Length = lengths,
                                Function  = annots_prot$mol_func
                                )

        left_annot <- rowAnnotation(df = annot_left,
                        col = list(Length = col_fun_len,
                                Function = color_prot_class),
                        simple_anno_size = unit(10, "mm"), 
                        width = NULL,
                        show_annotation_name=TRUE)
}



## TOP ANNOTATION
# SOURCE AND MOLECULAR FUNCTION
annot_top = data.frame(Source = annots_prot$source,
                       Function  = annots_prot$mol_func)

ha_top <- HeatmapAnnotation(df = annot_top,
            col = list(Source = color_source,
                        Function = color_prot_class),
            show_legend = c(TRUE, FALSE),
            simple_anno_size = unit(10, "mm"), width = NULL)



#file_fig_heat <- file.path(folder_path, 'heatmap_prots.png')

file_fig_heat <- paste0('test_heatmap_',database_name, '.png')

png(file = file_fig_heat, width = 1200, height = 1000)
Heatmap(mat_rmsd_proteins, 
        col = col_fun, name = 'Protein RMSD',
        clustering_distance_rows = "euclidean",
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        heatmap_legend_param = list(at = list_legend,
                                legend_direction = "vertical",
                                legend_height = unit(10, "cm") 
                                ),
        border = TRUE,
        top_annot = ha_top,
        left_annotation = left_annot,
        )
dev.off()

print('xxxxx')
print(database_name)
