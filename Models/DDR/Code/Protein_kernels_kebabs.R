## instead of user provided sequences in XStringSet format
## for this example a set of DNA sequences is created
## RNA- or AA-sequences can be used as well with the mismatch kernel

#Load Kebab
#library(kebabs, lib.loc = '/home/sevastopol/R/x86_64-pc-linux-gnu-library/3.6')
library(kebabs)
args = commandArgs(trailingOnly = TRUE)

generate_protein_kernels_using_KeBABs <- function(path_to_fasta,savepath){

    message('Reading file...')
    AAseq_fasta <- seqinr::read.fasta(path_to_fasta)

    format_AAseq_fasta <- function(AAseq_f){

        AAseq_f <- sapply(AAseq_f,function(x){ return(paste0(toupper(as.character(x)),collapse=''))})
        AAseq_c <- c()
        for (i in seq_along(AAseq_f)){
            AAseq_c <- c(AAseq_c,AAseq_f[[i]])
        }
        return(AAseq_c)
    }

    retrieve_AAsequences <- function(AAseq,namesSeq){
        dnaseqs <- AAStringSet(AAseq)
        #names(dnaseqs) <- paste("S", 1:length(dnaseqs), sep="")
        names(dnaseqs) <- namesSeq
        return(dnaseqs)
    }

    mismatch_kernel <- function(k,m,AAseq,namesSeq){
        mK <- mismatchKernel(k=k, m=m)
        return (mK(retrieve_AAsequences(AAseq,namesSeq)))
    }

    spectrum_kernel <- function(k,AAseq,namesSeq){
        mK <- spectrumKernel(k=k)
        return (mK(retrieve_AAsequences(AAseq,namesSeq)))
    }

    get_mismatch_kernels <- function(AAseq_fasta){

        if (!file.exists(paste0(savepath,'mismatch_kernel_k3m1.tsv'))){
            message('Computing mismatch kernel with k=3 and m=1')
            #k=3,m=1
            k3m1_Kernel <- mismatch_kernel(3,1,format_AAseq_fasta(AAseq_fasta),names(AAseq_fasta))
            message('Saving it!')
            write.table(k3m1_Kernel,paste0(savepath,'mismatch_kernel_k3m1.tsv'),sep='\t')
        }

        if (!file.exists(paste0(savepath,'mismatch_kernel_k4m1.tsv'))){
            message('Computing mismatch kernel with k=4 and m=1')
            #k=4,m=1
            k4m1_Kernel <- mismatch_kernel(4,1,format_AAseq_fasta(AAseq_fasta),names(AAseq_fasta))
            message('Saving it!')
            write.table(k4m1_Kernel,paste0(savepath,'mismatch_kernel_k4m1.tsv'),sep='\t')
        }

        if (!file.exists(paste0(savepath,'mismatch_kernel_k3m2.tsv'))){
            message('Computing mismatch kernel with k=3 and m=2')
            #k=3,m=2
            k3m2_Kernel <- mismatch_kernel(3,2,format_AAseq_fasta(AAseq_fasta),names(AAseq_fasta))
            message('Saving it!')
            write.table(k3m2_Kernel,paste0(savepath,'mismatch_kernel_k3m2.tsv'),sep='\t')
        }

        if (!file.exists(paste0(savepath,'mismatch_kernel_k4m2.tsv'))){
            message('Computing mismatch kernel with k=4 and m=2')
            #k=4,m=2
            k4m2_Kernel <- mismatch_kernel(4,2,format_AAseq_fasta(AAseq_fasta),names(AAseq_fasta))
            message('Saving it!')
            write.table(k4m2_Kernel,paste0(savepath,'mismatch_kernel_k4m2.tsv'),sep='\t')
        }

    }

    get_spectrum_kernel <- function(AAseq_fasta){

        if (!file.exists(paste0(savepath,'spectrum_kernel_k3.tsv'))){
            message('Computing spectrum kernel with k=3')
            #k=3
            k3_Kernel <- spectrum_kernel(3,format_AAseq_fasta(AAseq_fasta),names(AAseq_fasta))
            message('Saving it!')
            write.table(k3_Kernel,paste0(savepath,'spectrum_kernel_k3.tsv'),sep='\t')
        }

        if (!file.exists(paste0(savepath,'spectrum_kernel_k4.tsv'))){
            message('Computing spectrum kernel with k=4')
            #k=4
            k4_Kernel <- spectrum_kernel(4,format_AAseq_fasta(AAseq_fasta),names(AAseq_fasta))
            message('Saving it!')
            write.table(k4_Kernel,paste0(savepath,'spectrum_kernel_k4.tsv'),sep='\t')
        }

    }

    message('Generating mismatch kernels')
    get_mismatch_kernels(AAseq_fasta)

    message('Generating spectrum kernels')
    get_spectrum_kernel(AAseq_fasta)

}

#Generate Protein Kernels using KeBABS
message('Generating protein kernels using KeBABS')
generate_protein_kernels_using_KeBABs(args[1],args[2])