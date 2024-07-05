#' A function to adjust for near-site clusters for single-end data
#'
#' @param anno_path The path to the genome annotation file, which should be in
#' the .gtf format.
#' @param clusters A data.frame containing the near-site clusters regions to be
#' adjusted. The required columns are seqnames, start, end, strand.
#' @param mappedTSS A data.frame containing the mapped TSSs based on the
#' corresponding single-end data. The required columns are seqnames, tss,
#' strand, counts.
#' @param learning_max_dist An integer to control the adjustment learning
#' process. Any mapped TSSs which are outside of this range (bp) in the
#' downstream of an annotated TSS will be ignored in the adjustment learning.
#' The default is 600(bp).
#' @param learning_min_dist An integer to control the adjustment learning
#' process. Any mapped TSSs which are within this range (bp) in the
#' downstream of an annotated TSS will be ignored in the adjustment learning.
#' The default is 0(bp).
#' @param learning_methods A character to specify the method used to obtain
#' the adjustment distance. The input has to be either "weighted_mean" or
#' "median". The default is set to "weighted_mean".
#' @param gene_of_interest A character or a vector of characters, specifying
#' the genes that are of interest in terms of predicting TSS clusters.
#' @param gene_extension An integer specifying the length (bp) of the extended region
#' before the annotated 5' of a gene when predicting TSS clusters.
#' Defaults to 500 (bp).
#' @param ncore An integer specifying the number of cores used for computation.
#' Defaults to 1.
#'
#' @return \code{adjustTSS} outputs a list containing two elements.
#' 1) learnt_adj_dist is an integer specifying the learnt adjustment distances;
#' 2) adj_TSS_clusters is a data.frame for the adjusted TSS clusters. It has
#' these columns: seqnames, tss5prime, width, strand, org_start, org_end.
#' tss5prime means the location for the 5' end of the adjusted TSS clusters,
#' and width specifies their lengths. org_start and org_end specify the
#' corresponding genomic regions before adjustment.
#' @export adjustTSS
#'
#' @import parallel
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom IRanges IRanges
#' @importFrom stats median
#' @importFrom dplyr %>% distinct group_by summarize n filter select

adjustTSS <- function(anno_path,
                      clusters,
                      mappedTSS,
                      learning_max_dist = 600,
                      learning_min_dist = 0,
                      learning_methods = "weighted_mean",
                      gene_of_interest = NULL,
                      gene_extension = 500,
                      ncore = 1){

  # loading reference
  txdb <- makeTxDbFromGFF(anno_path)
  gene_list <- genes(txdb)
  #tx_list_by_gene <- transcriptsBy(txdb,by=c("gene"))
  exons_list_by_gene <- reduce(exonsBy(txdb, by="gene"))

  # TSS adjusting learning
  print("Start learning the adjustment distance...")
  mappedTSS <- mappedTSS[,c("seqnames", "tss",
                            "strand", "counts")]
  TSS_adj_learnt_all <- TSS_adj_learning(txdb = txdb,
                                         clusters = clusters,
                                         mappedTSS = mappedTSS,
                                         learning_max_dist = learning_max_dist,
                                         learning_min_dist = learning_min_dist,
                                         learning_methods = learning_methods,
                                         ncore = ncore)
  print(paste0("The learnt adjustment distance is ",TSS_adj_learnt_all," bp."))
  mappedTSS <- NULL

  # reduce gene size and cluster size
  clusters_gr <- makeGRangesFromDataFrame(df = clusters,
                                          keep.extra.columns = TRUE)

  ## 1) reduce to genes that are of interests
  if(is.null(gene_of_interest)){
    all_genes <- names(gene_list)
  }else{
    if(sum(!gene_of_interest%in%names(gene_list))==0){
      all_genes <- gene_of_interest
    }else{
      stop(" \"gene_of_interest\" is not valide. Please provide real gene names.")
    }
  }

  genes_tobe_adj <- gene_list[all_genes]

  ## 2) reduce both the clusters and genes to they have overlap
  resized_genes_tobe_adj <- resize(genes_tobe_adj,width(gene_list)+gene_extension,fix = "end")
  overlaps <- findOverlaps(resized_genes_tobe_adj, clusters_gr)
  #overlaps <- findOverlaps(genes_tobe_adj, clusters_gr)
  overlaps <- data.frame(overlaps)
  genes_idx <- unique(overlaps$queryHits)
  cluster_idx <- unique(overlaps$subjectHits)

  genes_tobe_adj <- genes_tobe_adj[genes_idx]
  clusters <- clusters[cluster_idx,]

  gene_list <- genes_tobe_adj
  ## out reduced results: reduced genes
  all_genes <- names(gene_list)
  ## out reduced results: TSS clusters_before_adj
  clusters_gr <- makeGRangesFromDataFrame(df = clusters,
                                          keep.extra.columns = TRUE)

  # TSS cluster adjusting
  print("Start adjusting near-site clusters regions...")

  adj_TSS_clusters <- parallel::mclapply(1:length(gene_list), function(iii){
    tryCatch({
      this_gene_region <- gene_list[iii]
      extended_this_gene_region <- resize(this_gene_region,width(this_gene_region) + gene_extension,fix = "end")
      this_gene_xn <- exons_list_by_gene[names(this_gene_region)]
      this_gene_clusters <- query_based_on_single_regionInput(single_region = extended_this_gene_region,
                                                              query_regions = clusters_gr)

      outs <- TSS_adjust_perGene(this_gene_region = this_gene_region,
                                 this_gene_clusters = this_gene_clusters,
                                 this_gene_xn = this_gene_xn,
                                 this_sample_leart_dist = TSS_adj_learnt_all,
                                 gene_extension = gene_extension)


      if (iii%%500 == 0){
        gc()
      }
      return(outs)
    },
    error=function(e){ })
  },
  mc.cores = ncore)

  print("Done!")

  adj_TSS_clusters <- adj_TSS_clusters[lengths(adj_TSS_clusters)>0]
  adj_TSS_clusters <- do.call(rbind,adj_TSS_clusters)

  return(list(learnt_adj_dist = TSS_adj_learnt_all,
              adj_TSS_clusters = adj_TSS_clusters)
         )
}
