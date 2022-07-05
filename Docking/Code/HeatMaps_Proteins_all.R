
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

## Load files
FILE_ANNOT <- "annot_test_chemlb.csv"
FILE_SIM <- "rmsd_test_chembl.csv"

# tests

# PROTEINS
# Mat rmsd protein
mat_rmsd_proteins <- read.table(FILE_SIM, sep = ";", header=T)
mat_rmsd_proteins <- as.matrix((mat_rmsd_proteins))
dim(mat_rmsd_proteins)
# Protein annotation
annots_prot <- read.table(FILE_ANNOT, header = TRUE, sep = ";", quote = "")




# color gradient of matix
col_fun <- colorRamp2(c(max(mat_rmsd_proteins), 20,5,0), c("#000000",  "#2a63a0", "#87add6","white"))
col_fun(seq(-100, 100))


# Load Paletes
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vec_super = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))


# Source
tags_source = unique(annots_prot$Source)
if (length(tags_source)==1) { 
    color_vec_source = c('#b80050')
    } else if (length(tags_source)==2) {
        color_vec_source = c('#b80050','#23074b')
        }else if (length(tags_source)>2) {
        color_vec_source = col_vec_super[1:length(tags_source)]
}

color_source = setNames(color_vec_source, tags_source)

### Protein class

tag_prot_class = unique(annots_prot$Class)[unique(annots_prot$Class) != "-"]


if (length(tag_prot_class)==1) { 
    color_vec_protclass = c('#b3e065')
    } else if (length(tag_prot_class)==2) {
        color_vec_protclass = c('#b3e065','#9a75cd')
        }else if (length(tag_prot_class)>2) {
        color_vec_protclass = col_vec_super[1+6:(length(tag_prot_class)+5)]
}


color_prot_class= setNames(color_vec_protclass, tag_prot_class)

color_prot_class <- c(color_prot_class, setNames('#ffffff', "-"))




df_anot_top = data.frame(source = annots_prot$Source[1:2000])

top_annot <- HeatmapAnnotation(df = df_anot_top,
        col = list(source = color_source),
        simple_anno_size = unit(10, "mm"), width = NULL)


bottom_annotation <- None
#
df_left_top = data.frame(mol_funct = annots_prot$Class[1:2000])

left_annot <- rowAnnotation(df = df_left_top,
                    col = list(mol_funct = color_prot_class),
        simple_anno_size = unit(10, "mm"), width = NULL)



pdf(file = 'test_chembl_heatmap_2.pdf', width=20, height=15)
Heatmap(mat_rmsd_proteins[1:2000,1:2000], col = col_fun, name = 'Protein RMSD',
        show_row_names = FALSE, show_column_names = FALSE,
        show_row_dend = FALSE, show_column_dend = FALSE,
        top_annotation = top_annot, 
        left_annotation=left_annot,
        border = TRUE
        )
dev.off()