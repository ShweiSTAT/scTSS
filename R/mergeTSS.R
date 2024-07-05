#' Merge TSS clusters for multiple samples
#'
#' @param unmerged_TSS_cluster_list A list of data.frames. Each data.frame
#' contains the TSS clusters for one sample. For paired-end data, each
#' data.frame has these columns: seqnames, start, end and strand. For TSS
#' clusters predicted from singe-end data, each data.frame has these columns:
#' seqnames, tss5prime, width, strand, org_start and org_end.
#' @param merge_type A character specifying the method used to merge
#' TSS clusters. There are two methods: "disjoint" and "reduce". Please refer
#' to \link[GenomicRanges:disjoin]{disjoin},
#' and \link[GenomicRanges:reduce]{reduce}
#' for more information. As for single-end data, only "disjoint"
#' is supported. Defaults to "disjoint".
#' @param if_paired FALSE/TRUE, specifying if we have paired-end data
#' (TRUE) or single-end data (FALSE). Default set to TRUE.
#' @param ncore An integer specifying the number of cores used for computation.
#' This functionality is disabled for the paired-end data. Defaults to 1.
#'
#' @return For paired-end data, \code{mergeTSS} outputs a data.frame, which
#' stores the unified TSS clusters across samples. The data.frame contains four
#' columns: seqnames, start, end and strand. For single-end data, the output is
#' a list, which contains two objects. 1) A data.frame containing the unified TSS
#' clusters across samples (columns: seqnames, tss5prime,
#' width and strand); 2) A list of data.frames. Each of them stores the
#' relations between unified TSS clusters and their corresponding
#' near-end genomic regions for a sample.
#' @export mergeTSS
#'
#' @import parallel
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomicFeatures


mergeTSS <- function(unmerged_TSS_cluster_list,
                     merge_type = "disjoin",
                     if_paired = TRUE,
                     ncore = 1){
  if(if_paired){
    unmerged_TSS_cluster_table <- do.call(rbind,unmerged_TSS_cluster_list)
    unmerged_TSS_cluster_list <- NULL

    unmerged_TSS_cluster_table <- makeGRangesFromDataFrame(df = unmerged_TSS_cluster_table,
                                                           keep.extra.columns = FALSE)
    if(merge_type == "disjoin"){
      merged_TSS_cluster_table <- disjoin(unmerged_TSS_cluster_table)
    }else{
      if(merge_type == "reduce"){
        merged_TSS_cluster_table <- reduce(unmerged_TSS_cluster_table)
      }else{
        stop("The merge_type has to be either \"reduce\" or \"disjoin\".")
      }
    }
    merged_TSS_cluster_table <- as.data.frame(merged_TSS_cluster_table)
    colnames(merged_TSS_cluster_table) <- c("seqnames",
                                            "start",
                                            "end",
                                            "width",
                                            "strand")
    merged_TSS_cluster_table <- merged_TSS_cluster_table[,c("seqnames",
                                                            "start",
                                                            "end",
                                                            "strand")]

    return(merged_TSS_cluster_table)

  }else{
    ## make the standard TSS cluster regions across samples
    unmerged_TSS_cluster_table <- do.call(rbind,unmerged_TSS_cluster_list)

    ## convert TSS clusters regions info to suit IRanges

    unmerged_TSS_cluster_table_sorted <- unmerged_TSS_cluster_table[,c("seqnames",
                                                                       "tss5prime",
                                                                       "width",
                                                                       "strand")]

    unmerged_TSS_cluster_table_sorted_plus <- unmerged_TSS_cluster_table_sorted[which(unmerged_TSS_cluster_table_sorted$strand=="+"),]
    unmerged_TSS_cluster_table_sorted_plus$start <-  unmerged_TSS_cluster_table_sorted_plus$tss5prime
    unmerged_TSS_cluster_table_sorted_plus$end <- unmerged_TSS_cluster_table_sorted_plus$tss5prime+unmerged_TSS_cluster_table_sorted_plus$width - 1

    unmerged_TSS_cluster_table_sorted_minus <- unmerged_TSS_cluster_table_sorted[which(unmerged_TSS_cluster_table_sorted$strand=="-"),]
    unmerged_TSS_cluster_table_sorted_minus$start <- unmerged_TSS_cluster_table_sorted_minus$tss5prime - unmerged_TSS_cluster_table_sorted_minus$width + 1
    unmerged_TSS_cluster_table_sorted_minus$end <- unmerged_TSS_cluster_table_sorted_minus$tss5prime


    unmerged_TSS_cluster_table_sorted <- rbind(unmerged_TSS_cluster_table_sorted_plus,
                                               unmerged_TSS_cluster_table_sorted_minus)
    unmerged_TSS_cluster_table_sorted <- makeGRangesFromDataFrame(df = unmerged_TSS_cluster_table_sorted,
                                                                  keep.extra.columns = FALSE)

    ## merge
    if(merge_type == "disjoin"){
      merged_TSS_cluster_table_sorted <- disjoin(unmerged_TSS_cluster_table_sorted)
    }else{
        stop("ERROR: this is single-end data. The merge_type has
             to be \"disjoin\".")
    }


    merged_TSS_cluster_table_out <- data.frame(seqnames = as.character(seqnames(merged_TSS_cluster_table_sorted)),
                                               tss5prime = start(promoters(merged_TSS_cluster_table_sorted, upstream = 0, downstream = 1)),
                                               width = width(merged_TSS_cluster_table_sorted),
                                               strand = as.character(strand(merged_TSS_cluster_table_sorted)))

    ## match the original near-site regions and merged TSS clusters
    merged_TSS_cluster_table_sorted <- as.data.frame(merged_TSS_cluster_table_sorted)
    merged_TSS_cluster_table_sorted <- cbind(merged_TSS_cluster_table_sorted,
                                             data.frame(tss5prime = merged_TSS_cluster_table_out$tss5prime))


    merged_TSS_cluster_table_sorted_gr <- makeGRangesFromDataFrame(df = merged_TSS_cluster_table_sorted,
                                                                   keep.extra.columns = TRUE)
    merged_TSS_cluster_table_per_sampe_out <- list()
    for( ii in 1:length(unmerged_TSS_cluster_list)){
      temp_unmerged_TSS_cluster_table <- unmerged_TSS_cluster_list[[ii]]

      temp_plus <- temp_unmerged_TSS_cluster_table[which(temp_unmerged_TSS_cluster_table$strand == "+"),]
      temp_plus$start <- temp_plus$tss5prime
      temp_plus$end <- temp_plus$tss5prime + temp_plus$width -1

      temp_minus <- temp_unmerged_TSS_cluster_table[which(temp_unmerged_TSS_cluster_table$strand == "-"),]
      temp_minus$start <- temp_minus$tss5prime - temp_minus$width +1
      temp_minus$end <- temp_minus$tss5prime

      temp_unmerged_TSS_cluster_table <- rbind(temp_plus,
                                               temp_minus)

      temp_unmerged_TSS_cluster_table_gr <- makeGRangesFromDataFrame(df = temp_unmerged_TSS_cluster_table,
                                                                     keep.extra.columns = TRUE)


      temp_rst <- parallel::mclapply(1:length(temp_unmerged_TSS_cluster_table_gr),
                                     function(iii){
                                       tryCatch({
                                       merge_match_for_SingleEnd_data(single_unmerged_region  = temp_unmerged_TSS_cluster_table_gr[iii],
                                                                      merged_regions = merged_TSS_cluster_table_sorted_gr)
                                       }, error=function(e){})
                                     },
                                     mc.cores = ncore)
      temp_rst <- do.call(rbind,temp_rst)
      merged_TSS_cluster_table_per_sampe_out[[ii]] <- temp_rst
      temp_rst <- NULL
    }
    names(merged_TSS_cluster_table_per_sampe_out) <- names(unmerged_TSS_cluster_list)

    return(list(merged_TSS_cluster_table = merged_TSS_cluster_table_out,
                cluster_to_region = merged_TSS_cluster_table_per_sampe_out))

  }
}
