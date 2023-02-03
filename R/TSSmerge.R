#' TSS merging for TSS count matrices of multiple samples
#'
#' @param tss_matrices A list of TSS count matrices output
#' from \code{TSSpredict}. Each element should a matrix.
#' @param tss_filtering_distance An integer specifying the threshold distance
#' for TSS candidates that should be considered from the same genuine TSS.
#' This should be set to the same as in \code{TSSpredict}. Defaults to 1.
#' @param ncore An integer specifying the number of cores used. Defaults to 1.
#'
#' @return \code{TSSmerge} outputs a matrix of read counts for the
#' predicted TSS by cell. This matrix is stored in data.table format.
#' @export TSSmerge
#' @importFrom  data.table .SD := data.table rbindlist set setDT
#' @import parallel
#'
#' @author Shiwei Fu, \email{shiwei.fu@email.ucr.edu}

TSSmerge <- function(tss_matrices,
                      tss_filtering_distance,
                      ncore = 1){

  ## step1: read in all samples and summarize all genes to be handled
    allgene <- character()
    gene_bentch <- list()
    for(i in 1:length(tss_matrices)){
      handle_outs <- tss_matrices[[i]]
      this_sample_gene <- unique(do.call(rbind, strsplit(handle_outs$unmerged_tss, ":", fixed=TRUE))[,2])
      gene_bentch[[i]] <- do.call(rbind, strsplit(handle_outs$unmerged_tss, ":", fixed=TRUE))[,2]
      allgene <- unique(c(allgene,this_sample_gene))
      handle_outs <- NULL
    }

    ## step2: extract unmerged TSS by gene
    unmerged_tss <- list()
    for( i in 1:length(allgene)){
      gene_id <- allgene[i]
      tss_names <- character()
      for(x in 1:length(tss_matrices)){
        idx <- which(gene_bentch[[x]] == gene_id)
        if(length(idx)>0){
          tss_names <- c(tss_matrices[[x]]$unmerged_tss[idx],tss_names)
        }
      }
      unmerged_tss[[i]] <- unique(tss_names)
    }

    ## step3: building merging dictionary
    merging_dictionary <- mclapply(unmerged_tss,function(tss){
      tryCatch({
        outs <- building_merging_dictionary(tss_names = tss,
                                            tss_filtering_distance = tss_filtering_distance)
        return(outs)
      },
      error=function(e){cat(" ERROR :",conditionMessage(e), "\n")})
    },mc.cores = ncore)
    merging_dictionary <- rbindlist(merging_dictionary,use.names=T)
    setDT(merging_dictionary,key="unmerged_tss")

    ## step4: update TSS information by dictionary
      for ( i in 1:length(tss_matrices)){
        single_sample <- tss_matrices[[i]]
        single_sample <- updateMatrix(mat = single_sample,
                                      merging_dictionary = merging_dictionary)
        tss_matrices[[i]] <- single_sample
        single_sample <- NULL
      }

    return(tss_matrices)


}
