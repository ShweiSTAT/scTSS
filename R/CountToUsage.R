#' Convert TSS count to TSS usage
#'
#' @param tss_count_matrix A TSS count matrix.
#' @param num_cell An integer specifying the threshold used to filter out
#' TSSs that are lowly expressed. A TSS must be expressed in this
#' many of cells for a single sample. Otherwise it will be removed.
#' Defaults to 30.
#'
#' @return \code{CountToUsage} outputs a matrix of usage for the
#' merged TSS by cell. This matrix is stored in data.table format.
#' @export CountToUsage
#' @importFrom data.table .SD := data.table rbindlist set setDT
#' @import parallel
#'
#' @author Shiwei Fu, \email{shiwei.fu@email.ucr.edu}

## ncore the num of core to run the methods
CountToUsage <- function(tss_count_matrix, num_cell=30){
  ##############
  ## version 2
  ##############
  ## remove lowly expressed TSS
  keep_idx <- which(rowSums(tss_count_matrix[,-1]>0) > num_cell)
  tss_count_matrix <- tss_count_matrix[keep_idx]

  ## remove genes that have only one TSS
  Gene_bench <- data.frame(do.call(rbind, strsplit(tss_count_matrix$merged_tss, ":", fixed=TRUE)))[,2]
  keep_idx <- which(Gene_bench %in% names(which(table(Gene_bench)>1)))
  tss_count_matrix <- tss_count_matrix[keep_idx]

  ## convert to usage
  Gene_bench <- data.frame(do.call(rbind, strsplit(tss_count_matrix$merged_tss, ":", fixed=TRUE)))[,2]
  tss_count_matrix <- cbind(Gene_bench,tss_count_matrix)

  tss_count_matrix <- tss_count_matrix[, (colnames(tss_count_matrix)[-c(1,2)]):=lapply(.SD, as.numeric), .SDcols = (colnames(tss_count_matrix)[-c(1,2)])]
  tss_usage_matrix <- tss_count_matrix[ ,(colnames(tss_count_matrix)[-c(1,2)]):=lapply(.SD, function (xx) xx/sum(xx)),
                                        .SDcols = (colnames(tss_count_matrix)[-c(1,2)]), by = "Gene_bench"]
  tss_count_matrix <- NULL
  removeNAs(tss_usage_matrix)
  tss_usage_matrix <- tss_usage_matrix[,-c("Gene_bench")]


  ##############
  ## version 1
  ##############
  # ## splice the the data set by gene
  # temp_inputSliced <- list()
  # Gene_bench <- data.frame(do.call(rbind, strsplit(tss_count_matrix$merged_tss, ":", fixed=TRUE)))[,2]
  # TobeGenes <- unique(Gene_bench)
  # for(ii in 1: length(TobeGenes)){
  #   TestGene_idx <- which(Gene_bench==TobeGenes[ii])
  #   temp_inputSliced[[TobeGenes[ii]]] <- tss_count_matrix[TestGene_idx]
  # }
  # tss_count_matrix <- NULL
  #
  # ## remove genes that have only one TSS
  # tss_num_of_gene <- unlist(lapply(temp_inputSliced, function(x){
  #   nrow(x)
  # }))
  # temp_inputSliced <- temp_inputSliced[which(tss_num_of_gene > 1)]
  #
  # ## remove lowly expressed Genes
  # expression_level_perCell <- mclapply(temp_inputSliced, function(x){
  #   colSums(x[,-1]>0)
  # },mc.cores = ncore)
  # gene_expression_level <- unlist(lapply(expression_level_perCell, function(x){
  #   sum(x>0)
  # }))
  # keep_idx <- which(gene_expression_level > num_cell)
  # expression_level_perCell <- expression_level_perCell[keep_idx]
  # temp_inputSliced <- temp_inputSliced[keep_idx]
  #
  # ## convert to usage
  # tss_usage_matrix <- mclapply(temp_inputSliced, function(x){
  #   ## sol-1
  #   x[ ,(colnames(x)[-1]):=lapply(.SD, function (xx) xx/sum(xx)),
  #            .SDcols = !"merged_tss"]
  #   removeNAs(x)
  #   return(x)
  # },mc.cores = ncore)
  # temp_inputSliced <- NULL
  # tss_usage_matrix <- rbindlist(tss_usage_matrix,use.names=T)

  return(tss_usage_matrix)
}
