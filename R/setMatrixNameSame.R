#' Prepare for the input of DU test
#'
#' @param Quantified_TSS_list A list of the outputs from the
#' "quantifyTSS" function. Each element of the list is a single output.
#' @param filtering_level A character, specifies how we filter out
#' lowly expressed TSS clusters for the DU test.
#' A TSS cluster is not necessarily
#' expressed in every sample. If we set the filtering_level to "cell" (default),
#' we only consider TSS clusters that are expressed
#' in more than a percentage (controlled by parameter `exp_level`)
#' of cells within each sample,
#' which is a conservative choice. If we set the
#' filtering_level to “sample”, we test TSS clusters that
#' are expressed in more than a certain percentage of samples.
#' This percentage is controlled by parameter `exp_level`.
#' This choice is more sensitive in finding DU TSS.
#' @param exp_level A percentage. When `filtering_level` is set to "cell,"
#' this parameter specifies the minimum percentage of cells within
#' a sample (as indicated by "sampleID" in the column metadata)
#' required for a TSS cluster to be considered valid.
#' When `filtering_level` is "sample",
#' it specifies the minimum percentage of samples
#' in which a TSS cluster must be expressed to be considered valid.
#' Default is 0, indicating we will not remove TSS clusters based on
#' expression levels.
#' @param remove_oneTSS_gene TRUE/FALSE, specifying if the genes with only one
#' TSS cluster will be removed. Default is TRUE, indicating the one-TSS
#' gene will be removed.
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
                              filtering_level="cell",
                              exp_level = 0,
                              remove_oneTSS_gene = TRUE
                              ){

  if(sum(!(filtering_level=="cell"|filtering_level=="sample"))>0){
    stop("Enter a correct filtering_level!" )
  }

  if(filtering_level=="cell"){
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
      print(paste0("In sample ", names(TSS_counts_list)[i]," ",nrow(temp) - length(temp_idx), " TSS clusters were removed for low expression."))
      TSS_counts_list[[i]] <- temp[temp_idx]
    }

    allnames <- data.table(TSS_clusters=Reduce(intersect,lapply(TSS_counts_list, function(y){y$TSS_clusters})))

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
                   " TSS clusters were removed for being the single TSS cluster on gene."))
    }

    ## outs
    genes_bentch <- do.call(rbind,
                            strsplit(TSS_counts_out$TSS_clusters,
                                     split = ":",
                                     fixed = TRUE))[,2]
    genes_retained <- unique(genes_bentch)
    print(paste0(nrow(TSS_counts_out), " TSS clusters on ",
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

  }else{
    TSS_counts_list <- list()
    col_meta_list <- list()
    for (i in 1:length(Quantified_TSS_list)) {
      temp <- Quantified_TSS_list[[i]]
      TSS_counts_list[[i]] <- temp$TSS_cluster_counts
      col_meta_list[[i]] <- temp$Column_meta
    }

    names(TSS_counts_list) <- names(Quantified_TSS_list)
    names(col_meta_list) <- names(Quantified_TSS_list)
    Quantified_TSS_list <- NULL

    ## step 1: Remove TSS clusters that have a low expression level
    # pre removal
    for (i in 1:length(TSS_counts_list)) {
      temp <- TSS_counts_list[[i]]
      temp_idx <- which(rowSums(temp[,-1] >0 ) > 1)
      #print(paste0("In sample ", names(TSS_counts_list)[i]," ",nrow(temp) - length(temp_idx), " TSS clusters were removed for low expression."))
      TSS_counts_list[[i]] <- temp[temp_idx]
    }
    # Actual removal: based on how often the TSS showed up samples
    TSS_cluster_freq <- table(unlist(lapply(TSS_counts_list, function(x){
      x$TSS_clusters
    })))

    TSS_retained <- names(which(TSS_cluster_freq/length(TSS_counts_list) > exp_level))

    allTSS <- unique(unlist(lapply(TSS_counts_list, function(y){y$TSS_clusters})))
    print(paste0("Across all samples, there are ", length(allTSS), " TSS clusters."))
    print(paste0(length(allTSS) - length(TSS_retained), " of them were removed for low expression."))
    allTSS <- TSS_retained

    # Actual removal: remove one-TSS genes
    if(remove_oneTSS_gene){
      genes_bentch <- do.call(rbind,
                              strsplit(allTSS,
                                       split = ":",
                                       fixed = TRUE))[,2]
      gene_freq <- table(genes_bentch)
      genes_retained <- names(which(gene_freq >1))
      allTSS <- allTSS[which(genes_bentch%in%genes_retained)]

      print(paste0(length(genes_bentch) - length(allTSS),
                   " TSS clusters were removed for being the single TSS cluster on gene."))
    }

    ## step 2: merge TSS clusters
    allnames <- data.table(TSS_clusters=allTSS)

    TSS_counts_out <- lapply(TSS_counts_list, function(y){
      y <- merge(x=allnames,y=y,by="TSS_clusters",all.x=TRUE)
      removeNAs(y)
      setDT(y,key="TSS_clusters")
    })

    TSS_counts_out <- Reduce(function(x, y) merge(x, y, by = "TSS_clusters", all = TRUE), TSS_counts_out)
    col_meta_out <- do.call(rbind,col_meta_list)
    rownames(col_meta_out) <- NULL

    ## outs
    genes_bentch <- do.call(rbind,
                            strsplit(TSS_counts_out$TSS_clusters,
                                     split = ":",
                                     fixed = TRUE))[,2]
    genes_retained <- unique(genes_bentch)
    print(paste0(nrow(TSS_counts_out), " TSS clusters on ",
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



}
