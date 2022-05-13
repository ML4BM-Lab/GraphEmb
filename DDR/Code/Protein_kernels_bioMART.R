#For this database, we went to https://m.ensembl.org/biomart/martview/
#and select all the GOs + gene names

#now we will process it through the R package csbl.go to obtain the adjacency matrix we are looking for
library('csbl.go')
args = commandArgs(trailingOnly = TRUE)

init_GO_list <- function(protpath){
    
    #init path
    bioMartDB <- paste0(getwd(),'/DBs/BioMART/mart_export.txt')

    #read file
    bioMARTf <- read.delim(bioMartDB)

    #select only desired gene-names
    selection_p <- list.files(protpath)

    message('Generating GO list')
    GO_list <- list()
    selectedfiles <- c()
    #iterate through all the files and save independently each file
    for (file in selection_p){
        #read symbols
        tryCatch(expr = {
            #read all symbols
            symbolf <- as.character(read.table(paste0(protpath,'/',file))[['V1']])
            #check which matches
            symbol_match <- names(which(sapply(symbolf,function(x) x%in%bioMARTf[['Gene.name']]))[1])
            if (!is.na(symbol_match)){
                #get the GO terms from the symbol that matches
                GOterms_file <- as.character(bioMARTf[which(bioMARTf[['Gene.name']]==symbol_match),'GO.term.accession'])
                #if its not empty
                if (length(GOterms_file)){
                    #save that as a GO list
                    GO_list[[substr(file,1,nchar(file)-4)]] <- paste0(c(symbol_match,GOterms_file),collapse=' ')
                    selectedfiles <- c(selectedfiles,file)
                }
            }
        },error = function(e){
            #this is in case symbolf returns an error becasue the file is empty
            message('Empty file')
        })
    }

    files_to_iterate <- sort(as.character(sapply(selectedfiles, function(x) substr(x,1,nchar(x)-4))))
    return(list(GOl=GO_list,fti=files_to_iterate))
}

##GUIDE OF CSBL.GO
#https://csbi.ltdk.helsinki.fi/csbl.go/userguide.pdf (MOLECULAR FUNCTION)

## CSBL.GO
get_sim_mat_from_csblgo <- function(dataset,savepath,GO_list,files_to_iterate){

    #save GO list item by item
    savep <- paste0('/tmp/bioMART_export_',dataset,'_prepmat.txt')

    if (file.exists(savep)) {
        #Delete file if it exists
        file.remove(savep)
    }
    #save
    for (hsagene in files_to_iterate){
        write(GO_list[[hsagene]],savep,append=TRUE)
    }

    #prepare matrix
    ent <- entities.from.text(savep)

    #delete those with MF = 0
    MFlen <- sapply(ent$MF,length)
    MFnonzero <- which(!MFlen == 0)
    ent <- ent[MFnonzero]
    #update files according to that condition
    files_to_iterate <- files_to_iterate[MFnonzero]

    #assign probabilities
    set.prob.table(organism=TAXONOMY.HUMAN, type="similarity")

    #compute similarity matrix
    message('Generating MF GO adj matrix of shape ',length(files_to_iterate), ' square matrix')
    sim_mat_PPI_MF <- entity.sim.many(ent, "MF", "Resnik")
    message('Done!')
    # sim_mat_PPI_BP <- entity.sim.many(ent, "BP", "Resnik")
    # sim_mat_PPI_CC <- entity.sim.many(ent, "CC", "Resnik")

    #check NA
    check_na <- table(is.na(sim_mat_PPI_MF))
    message('NANs Table: ',check_na['FALSE'], ' false and ',check_na['TRUE'], ' true')

    #normalize
    sim_mat_PPI_MF_norm <- sim_mat_PPI_MF/apply(sim_mat_PPI_MF,1,max)

    #colnames & rownames
    colnames(sim_mat_PPI_MF_norm) <- files_to_iterate
    rownames(sim_mat_PPI_MF_norm) <- files_to_iterate

    #save
    message('Saving..')
    write.table(sim_mat_PPI_MF_norm, paste0(savepath,'GO_PPI.tsv'),sep='\t')
}

#call and generate the PPI GO matrix
init_l <- init_GO_list(protpath=args[1])
GO_l <- init_l[['GOl']]
files_to_iterate <- init_l[['fti']]
message('Generating simmat using csblgo R library')
get_sim_mat_from_csblgo(dataset=args[2],savepath=args[3],GO_l,files_to_iterate)

