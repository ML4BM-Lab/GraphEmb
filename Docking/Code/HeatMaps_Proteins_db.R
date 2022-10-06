
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
database_name <- 'BIOSNAP'

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


##############
## PLOT FIRST HISTOGRAM
hits_prot = mat_rmsd_proteins[upper.tri(mat_rmsd_proteins, diag = FALSE)]
data_hist <- data.frame(hits_prot)

file_fig_hist <- paste0(folder_db_path, '/hist_prot_', database_name, '.pdf')

file_fig_hist <- paste0('0_test_hist_prot_', database_name, '.pdf')

pdf(file = file_fig_hist, width=16, height=8)
ggplot(data_hist, aes(hits_prot)) +
  ggtitle("Histogram of RMSD distribution")+
  labs(x = "RMSD (Å)") +             
  geom_histogram(bins=160, colour="#0c0f0a", fill="#f18f01", alpha=0.95, size=0.01)+
#  geom_density()+
  #xlim(0,40)+
  scale_x_continuous(limits = c(0,40),
                     breaks = seq(0,40,5))+
  scale_y_continuous(breaks=seq(0,100000,5000))+
  #scale_x_continuous(breaks=seq(0,40,5))+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size=20))   
dev.off()


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
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


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
tag_prot_class = unique(annots_prot$mol_func)

if (length(tag_prot_class)==1) { 
    color_vec_protclass = c('#b3e065')
    } else if (length(tag_prot_class)==2) {
        color_vec_protclass = c('#b3e065','#9a75cd')
        }else if (length(tag_prot_class)>2) {
        color_vec_protclass = col_vec_super[1+9:(length(tag_prot_class)+8)]
}

color_prot_class= setNames(color_vec_protclass, tag_prot_class)

#color_prot_class <- c(color_prot_class, setNames('#ffffff', "-"))

## Enzyme code
ec_code <- annots_prot$EC
tags_ec <- sort(unique(ec_code)[unique(ec_code) != "-"]) # change again for -
if (length(tags_ec)!=0){
    color_list_ec = c('#ea481f', '#ffad0a', '#6dcd3c', '#f50093', '#009eb4', '#6a1145', '#041f1e')
    color_ec <- setNames(color_list_ec[1:length(tags_ec)], tags_ec)
    color_ec <- c(color_ec, setNames('#ffffff', "-"))
}


#### Create annotations

## TOP ANNOTATION
# SOURCE AND MOLECULAR FUNCTION
annot_top = data.frame(Source = annots_prot$source,
                       Function  = annots_prot$mol_func)

ha_top <- HeatmapAnnotation(df = annot_top,
        col = list(Source = color_source,
                    Function = color_prot_class),
        show_legend = c(TRUE, FALSE),
        simple_anno_size = unit(10, "mm"), width = NULL)


# LENGTH AND FUNCTION
# min(annots_prot$length) changed for 0
# max(annots_prot$length) changed for 1000
lengths <- as.double(annots_prot$length)

col_fun_len <- colorRamp2(c(max(annots_prot$length), 800, 400,200, min(annots_prot$length)), 
c('#3d1e2e','#5d2e46', '#bc928a','#d0ada7', '#eedfd6'))
col_fun_len(seq(-150, 150))


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

### SAVE

#file_fig_heat <- file.path(folder_path, 'heatmap_prots.png')

file_fig_heat <- paste0(folder_db_path, '/heatmap_prots_', database_name, '.png')

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
