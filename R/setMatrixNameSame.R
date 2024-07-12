#' Prepare for the input of DU test
#'
#' @param Quantified_TSS_list A list of the outputs from the
#' "quantifyTSS" function. Each element of the list is a single output.
#' @param exp_level A percentage specifying, within a sample
#' (specified by "sampleID" in the column meta data) the minim percentage of
#' cells that should have a valid TSS cluster expressed. Default is 0,
#' indicating we will not remove TSS clusters based on expression levels.
#' @param remove_oneTSS_gene TRUE/FALSE, specifying if the genes with only one
#' TSS cluster will be removed. Default is TRUE, indicating the one-TSS
#' gene will be removed.
#' @param join_type a character, specifies how we prepare the TSS clusters
#' for the DU test. A TSS cluster is not necessarily
#' expressed in every sample. If we set the join_type to "inner",
#' we only consider TSS clusters that are expressed in all of the
#' samples, which is a conservative choice. If we set the
#' join_type to “outer”, we test TSS regions that
#' are expressed in as few as only one sample. This choice is more
#' sensitive in finding DU TSS. The default is "inner".
#'
#' @return \code{setMatrixNameSame} returns a list containing these two elements:
#' 1) TSS_count_matrix, a data.matrix specifying the TSS cluster counts
#' ready for DU test; 2) col_meta, a column meta data.frame compatible with
#' TSS_count_matrix, which is ready for DU test.
#' @export SetMatrixNameSame
#'
#' @importFrom data.table  .SD := data.table rbindlist set setDT
#'
SetMatrixNameSame <- function(Quantified_TSS_list,
                              exp_level = 0,
                              remove_oneTSS_gene = TRUE,
                              join_type="inner"){

  TSS_counts_list <- list()
  col_meta_list <- list()
  for (i in 1:length(Quantified_TSS_list)) {
    temp <- Quantified_TSS_list[[i]]
    TSS_counts_list[[i]] <- temp$TSS_cluster_counts
    col_meta_list[[i]] <- temp$Column_meta
  }

  names(TSS_counts_list) <- names(Quantified_TSS_list)
  names(col_meta_list) <- names(Quantified_TSS_list)

  ## step 1: remove TSS clusters that have a low expression level
  for (i in 1:length(TSS_counts_list)) {
    temp <- TSS_counts_list[[i]]
    temp_idx <- which(rowSums(temp[,-1] >0 ) / (ncol(temp) - 1 ) > exp_level)
    print(paste0("In sample ", names(TSS_counts_list)[i]," ",nrow(temp) - length(temp_idx), " TSS clusters were delected for low expression."))
    TSS_counts_list[[i]] <- temp[temp_idx]
  }

  ## step 2: merge TSS clusters
  if(sum(!(join_type=="inner"|join_type=="outer"))>0){
    stop("Enter a correct joint_type." )
  }

  if(join_type == "inner"){
    allnames <- data.table(TSS_clusters=Reduce(intersect,lapply(TSS_counts_list, function(y){y$TSS_clusters})))
  }else{
    allnames <- data.table(TSS_clusters=unique(unlist(lapply(TSS_counts_list, function(y){y$TSS_clusters}))))
  }

  TSS_counts_out <- lapply(TSS_counts_list, function(y){
    y <- merge(x=allnames,y=y,by="TSS_clusters",all.x=TRUE)
    removeNAs(y)
    setDT(y,key="TSS_clusters")
  })

  TSS_counts_out <- Reduce(function(x, y) merge(x, y, by = "TSS_clusters", all = TRUE), TSS_counts_out)
  col_meta_out <- do.call(rbind,col_meta_list)
  rownames(col_meta_out) <- NULL

  ## step 3: remove one-TSS gene
  if(remove_oneTSS_gene){
    genes_bentch <- do.call(rbind,
                            strsplit(TSS_counts_out$TSS_clusters,
                                     split = ":",
                                     fixed = TRUE))[,2]
    gene_freq <- table(genes_bentch)
    genes_retained <- names(which(gene_freq >1))
    TSS_counts_out <- TSS_counts_out[which(genes_bentch%in%genes_retained)]
    print(paste0(length(genes_bentch) - nrow(TSS_counts_out),
                 " TSS clusters were delected for being the single TSS cluster on gene."))
  }

  ## outs
  genes_bentch <- do.call(rbind,
                          strsplit(TSS_counts_out$TSS_clusters,
                                   split = ":",
                                   fixed = TRUE))[,2]
  genes_retained <- unique(genes_bentch)
  print(paste0(nrow(TSS_counts_out), " TSS clsuters on ",
               length(genes_retained), " genes are in the final output."))


  TSS_counts_out <- as.matrix(TSS_counts_out, rownames=1)
  TSS_counts_out <- TSS_counts_out[,col_meta_out$barcode]

  ## final check: make sure the columns are matched
  if(sum(colnames(TSS_counts_out) != col_meta_out$barcode) > 0){
    stop("Make sure the column names in the TSS clusters count matrices and
         the \'barcode\' column in the column meta data are matched")
  }

  return(list(TSS_count_matrix = TSS_counts_out,
              col_meta = col_meta_out))
}
