#' Make Mapped table for a sample
#'
#' @param anno_path The path to the genome annotation file, which should be in
#' the .gtf format.
#' @param bam_path The path to the .bam, which should be stored under the same
#' directory as its .bai file.
#' @param gene_of_interest A character or a vector of characters, specifying
#' the gene names that are of interest in terms of predicting Mapped TSS. The
#' input gene names should be same as those in the genome annotation
#' (see anno_path). The default is NULL, which means we need to check all genes
#' provided by the genome annotation (see anno_path).
#' @param if_paired FALSE/TRUE, specifying if we have paired-end (on-site) data
#' (TRUE) or single-end (near-site) data (FALSE). Default set to TRUE.
#' @param gene_extension An integer specifying the length of the extended region
#' before the annotated 5' of a gene when extracting mapped TSS. Defaults to 500 (bp).
#' @param ncore An integer specifying the number of cores used for computation.
#' Defaults to 1.
#'
#' @return \code{makeMappedTSS} outputs a data.frame of Mapped TSS table, which
#' contains four columns: seqnames, tss (the location of the Mapped TSS),
#' strand and counts (the count of reads aligned with the Mapped TSS)
#' @export makeMappedTSS
#'
#' @import parallel
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom Rsamtools ScanBamParam

makeMappedTSS <- function(anno_path,
                          bam_path,
                          gene_of_interest = NULL,
                          if_paired = TRUE,
                          gene_extension = 500,
                          ncore = 1){

  ## load in the needed files from path
  txdb <- makeTxDbFromGFF(anno_path)
  ## get info for genes to be tabulated
  gene_list <- genes(txdb)
  if(is.null(gene_of_interest)){
    all_genes <- names(gene_list)
  }else{
    if(sum(!gene_of_interest%in%names(gene_list))==0){
      all_genes <- gene_of_interest
    }else{
      stop("gene_of_interest is not valide. Please provide real gene names.")
    }
  }
  gene_list <- gene_list[all_genes]

  ## making Mapped TSS table for all the genes of interest

  print(paste0("Start making Mapped TSS table for ", length(gene_list)," genes") )

  MappedTSS_outs_all <- parallel::mclapply(1:length(gene_list), function(iii){
    tryCatch({
      gene_region_of_this_gene <- gene_list[iii]
      outs <- getMappedTSS_perGene(region_of_this_gene = gene_region_of_this_gene,
                                 bam_path = bam_path,
                                 if_paired_this_data = if_paired,
                                 gene_extension = gene_extension)
      if (iii%%500 == 0){
        gc()
      }
      return(outs)
    },
    error=function(e){cat("Gene: ",all_genes[iii]," ERROR :",conditionMessage(e), "\n")})
  },
  mc.cores = ncore)

  MappedTSS_outs_all <- MappedTSS_outs_all[which(lengths(MappedTSS_outs_all)>0)]
  MappedTSS_outs_all <- do.call(rbind,MappedTSS_outs_all)

  print("Finished making Mapped TSS table!")
  return(MappedTSS_outs_all)
}

