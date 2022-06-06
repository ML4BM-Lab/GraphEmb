## From the official source -> gklambauer/Rchemcpp

    # sd2gramSpectrum(sdf, sdf2,
    #          kernelType = c("spectrum", "tanimoto", "minmaxTanimoto", "marginalized", "lambda"),
    #          margKernelEndProbability = 0.1, lambdaKernelLambda = 1,
    #          depthMax = as.integer(3), onlyDepthMax = FALSE,
    #          flagRemoveH = FALSE, morganOrder = as.integer(0),
    #          silentMode = FALSE, returnNormalized = TRUE,
    #          detectArom = FALSE)

    # sdf: File containing the molecules. Must be in MDL file format
    #           (MOL and SDF files). For more information on the file format
    #           see http://en.wikipedia.org/wiki/Chemical_table_file.

#Load library (from other path)
library(Rchemcpp,lib.loc = '/home/sevastopol/R/x86_64-pc-linux-gnu-library/3.6')
library(BiocParallel)
args = commandArgs(trailingOnly = TRUE)

generate_drug_kernels_using_Rchemcpp <- function(path_to_molfiles, path_to_output, recompute=FALSE){
    
    #Load files
    sdf <- sort(list.files(path_to_molfiles,full.names=TRUE))

    #delete the empty ones
    for (sfile in sdf){
        SDFsetobject <- ChemmineR::read.SDFset(sfile)
        if (!ChemmineR::validSDF(SDFsetobject)){
            message(sfile, ' invalid')
            sdf <- setdiff(sdf,sfile)
        }
    }

    #now recomputes
    sdf_L <- length(sdf)

    #list of matrices
    sdf_matrix <- matrix(1,sdf_L,sdf_L)

    #total kernels
    tot_kernels <- c('spectrum','tanimoto','minmaxTanimoto','marginalized','lambda')

    #check if already exists
    for (kerneln in tot_kernels){
        if (file.exists(paste0(path_to_output,'Rchemcpp_',kerneln,'.tsv'))){
            
            #drop that kerneln if recompute=FALSE
            if (!recompute){
            message('kernel ',kerneln,' already exists, wont be computed')
            tot_kernels <- setdiff(tot_kernels,kerneln)
            #do not drop that kernel
            }else{
                message('kernel ',kerneln,' already exists, but as recompute is set to TRUE, it will be recomputed')
            }
        }
    }

    ComputeKernel <- function(ijelements,sdf,kerneln){
        i <- ijelements[1]
        j <- ijelements[2]
        return(sd2gramSpectrum(sdf[i],sdf[j],kernelType = kerneln,silentMode=TRUE))
    }

    message('Generating combinations...')
    #Generate combinations of files
    ijelementslol <- lapply(seq(sdf_L-1),function(i){
        lapply(seq(i+1,sdf_L),function(j){
                return (c(i,j))
        })
    })

    # ijelementslol <- unlist(lapply(seq(sdf_L-1),function(i){
    #     lapply(seq(i+1,sdf_L),function(j){
    #             return (c(i,j))
    #     })
    # }),recursive=FALSE)

    message('Computing kernels! Setting up the parallelization, using 60 cores')
    parallClass <- BiocParallel::bpparam()
    parallClass$workers <- 60

    #List of kernels
    for (kerneln in tot_kernels){
        message('Kernel -> ', kerneln)

        uppertriv <- c()
        row_i <- 0
        for (ijelementsl in ijelementslol){
            row_i <- row_i + 1
            message('Row number ',row_i)
            uppertriv <- c(uppertriv,unlist(BiocParallel::bplapply(ijelementsl,ComputeKernel,sdf,kerneln)))
        }
        # for (ijelementsl in ijelementslol){
        #     row_i <- row_i + 1
        #     message('Row number ',row_i)
        #     uppertriv <- c(uppertriv,ComputeKernel(ijelementsl[[1]],sdf,kerneln))
        # }

        pp <- sapply(ijelementsl,function(el){
            print(el)
            ComputeKernel(el,sdf,kerneln)
        })

        # uppertriv <- unlist(BiocParallel::bplapply(ijelementslol,ComputeKernel,sdf,kerneln))
        sdf_matrix[upper.tri(sdf_matrix)] <- uppertriv
        sdf_matrix[lower.tri(sdf_matrix)] <- sdf_matrix[upper.tri(sdf_matrix)]
        diag(sdf_matrix) <- rep(1,sdf_L)

        #Save file
        message('Saving ',kerneln)
        write.table(sdf_matrix,paste0(path_to_output,'Rchemcpp_',kerneln,'.tsv'),sep='\t')
    }
}


#Generate Drug Kernels using Rchemcpp
message('Generating Drug kernels using Rchemcpp, with args ',args)
generate_drug_kernels_using_Rchemcpp(args[1],args[2])
