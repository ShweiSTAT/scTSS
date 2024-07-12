#' DU TSS test
#'
#' @param mat A data.matrix for the count data of TSS clusters.
#' The rows are TSS clusters and the columns are cells.
#' @param col_meta A data.frame containing the column meta information for "mat".
#' It has to contain at least these columns: barcode, condition, sampleID.
#' @param agg_level A string to specify the DU test model ("bulk" or "cell").
#' "cell" means the DU test are conducted at single-cell level,
#' "bulk"means the DU test are conducted at the pseudo-bulk level.
#' The defult is "bulk."
#' @param full_model A formula for the full model. The default is
#' cbind(TSS_counts,gene_counts - TSS_counts) ~ condition + (1|sampleID).
#' @param base_model A formula for the base model. The default is
#' cbind(TSS_counts,gene_counts - TSS_counts) ~ (1|sampleID).
#' @param ncore An integer specifying the number of cores used for computation.
#' Defaults to 1.
#'
#' @return \code{DUtest} returns a data.frame, containing the test results for
#' each TSS cluster.
#' @export DUtest
#'
#' @import parallel
#' @import lme4
#' @importFrom stats coef p.adjust
#'
DUtest <- function(mat,
                   col_meta,
                   agg_level = "bulk",
                   full_model = cbind(TSS_counts,gene_counts - TSS_counts) ~ condition + (1|sampleID),
                   base_model = cbind(TSS_counts,gene_counts - TSS_counts) ~ (1|sampleID),
                   ncore = 1){

  ## pre-step 1: check if col_meta and the mat have the same name

  if(sum(col_meta$barcode != colnames(mat)) >0){
    stop("The columns in the TSS count matrix and the col_meta should be matched.")
  }

  ## pre-step 2: check if all the columns are included in the meta data
  if(sum(c("barcode","condition","sampleID" )%in%colnames(col_meta) ) !=3){
    stop("Please make sure \'barcode\',\'condition\',\'sampleID\' are included in all column meta data.")
  }


  ## DU test
  print(paste0("Start running DU test."))
  # slice TSS counts by genes
  rowname_bentch <- do.call(rbind,strsplit(rownames(mat),":",fixed = TRUE))[,2]
  gene_all <- unique(rowname_bentch)
  sliced_TSS_counts <- list()
  for (i in 1:length(gene_all)){
    idx <- which(rowname_bentch == gene_all[i])
    temp <- mat[idx,]
    sliced_TSS_counts[[gene_all[i]]] <- temp
  }
  mat <- NULL

  # Test by gene
  rst <- parallel::mclapply(1:length(sliced_TSS_counts), function(x){
    tryCatch({
      suppressWarnings(suppressMessages(outs <- glmm_test_for_a_gene(temp_full_data = sliced_TSS_counts[[x]],
                                   col_meta = col_meta,
                                   agg_level = "bulk",
                                   full_model = full_model,
                                   base_model = base_model)))

      return(outs)
    },
    error=function(e){ })

  },mc.cores = ncore)

  outs <- do.call(rbind,rst)
  # outs$pval_adj <- p.adjust(outs$pval,
  #                           method = "BH")

  print("Done!")
  return(outs)



}
