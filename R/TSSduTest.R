#' DU TSS test
#'
#' @param DU_input A list of TSS usage matrices.
#' The order the matrices should match the sample label.
#' @param sample_label A vector of characters specifying the condition
#' labels for the input matrices. The number of input matrices and vector
#' of labels should be of the same length.
#' @param permutation_times An integer specifying the permutation times
#' for DU test. Defaults to 2000.
#' @param join_type An character, either being "inner" or "outer",
#' specifying how we unify the merged TSS regions across input usage matrices.
#' "inner" means we only keep the TSS regions that are shared
#' by all input matrices for differential usage test.
#' "outer" means we keep all TSS regions for all input matrices.
#' @param ncore An integer specifying the number of cores used. Defaults to 1.
#'
#' @return \code{TSSduTest} outputs the DU test results and the
#' design matrix used for DU test.
#' @export TSSduTest
#' @import data.table
#' @import parallel
#' @importFrom  stats aggregate median model.matrix
#' @importFrom  transport wasserstein1d
#'
#' @author Shiwei Fu, \email{shiwei.fu@email.ucr.edu}

TSSduTest <- function(DU_input,
                      sample_label,
                      permutation_times=2000,
                      join_type = "inner",
                      ncore=1){
  ## the number of permutations
  a <- permutation_times
  ## sample label
  sample_label <- data.frame(sample_label=sample_label)
  #############################
  ### preprocessing before DU
  #############################
  ## set TSS names same across samples
  DU_input <- setTSSnamesSame(x=DU_input,
                              join_type = join_type)
  ## make the covariate matrix
  z <- model.matrix(~-1 + sample_label, data = sample_label)
  z <- as.matrix(z[,-1])
  ## permute covariate matrix
  z_per <- list()
  set.seed(1)
  for( i in 1:a){
    z_per[[i]] <- z[sample(nrow(z)),]
  }
  ## slice input usage by genes
  DU_inputSliced <- list()
  Gene_bench <- data.frame(do.call(rbind, strsplit(DU_input[[1]]$merged_tss, ":", fixed=TRUE)))[,2]
  TobeGenes <- names(which(table(Gene_bench) > 1))
  for(i in 1: length(TobeGenes)){
    TestGene_idx <- which(Gene_bench==TobeGenes[i])
    DU_inputSliced[[TobeGenes[i]]] <- lapply(DU_input, function(x){
      return(x[TestGene_idx])
    })
  }
  DU_input <- NULL

  ###########
  ## DU test
  ###########
  DU_outs <- parallel::mclapply(1:length(DU_inputSliced),function(iii){
    tryCatch({
      outs <- DU_test(test_gene = DU_inputSliced[[iii]],
                      a = a,
                      z = z,
                      z_per = z_per,
                      sample_label = sample_label)
      return(outs)
    },
    error=function(e){cat(iii,"Gene: ",TobeGenes[iii]," ERROR :",conditionMessage(e), "\n")})
  },mc.cores = ncore)
  DU_outs <- do.call(rbind,DU_outs)

  return(list(DU_result = DU_outs,
              design_matrix = z))

}
