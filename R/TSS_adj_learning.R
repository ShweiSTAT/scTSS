# a function to learn tss adjust distance for a sample

TSS_adj_learning <- function(txdb,
                             clusters,
                             mappedTSS,
                             learning_max_dist = 600,
                             learning_min_dist = 0,
                             learning_methods = "weighted_mean",
                             ncore = 1){

  if(! (learning_methods=="median"|learning_methods == "weighted_mean") ){
    stop("learning_methods is not valide. Please choose between 'median' or 'weighted_mean'.")
  }

  gene_list <- genes(txdb)
  exons_list_by_gene <- reduce(exonsBy(txdb, by="gene"))
  clusters_gr <- makeGRangesFromDataFrame(df = clusters,
                                          keep.extra.columns = TRUE)
  mappedTSS_gr <- GRanges(seqnames = mappedTSS$seqnames,
                          ranges = IRanges(start = mappedTSS$tss, end = mappedTSS$tss),
                          strand = mappedTSS$strand)
  ##################
  ## prep
  ##################

  ## one TSS gene name if we use one exon genes to learn
  oneTSSgene <- extract_oneTSS_gene_from_Anno(txdb)
  onTSSgene_exons <- exons_list_by_gene[oneTSSgene]
  genes_tobe_learnt <- names(which(lengths(onTSSgene_exons) == 1))

  ## reduce the size of data to genes only to be learnt
  genes_tobe_learnt <- gene_list[genes_tobe_learnt]
  # reduce the clusters and genes
  overlaps <- findOverlaps(genes_tobe_learnt, clusters_gr)
  overlaps <- data.frame(overlaps)
  genes_idx <- unique(overlaps$queryHits)
  cluster_idx <- unique(overlaps$subjectHits)

  genes_tobe_learnt <- genes_tobe_learnt[genes_idx]
  clusters <- clusters[cluster_idx,]

  # reduce the Mapped TSS
  overlaps <- findOverlaps(genes_tobe_learnt, mappedTSS_gr)
  overlaps <- data.frame(overlaps)
  genes_idx <- unique(overlaps$queryHits)
  mappedTSS_idx <- unique(overlaps$subjectHits)

  genes_tobe_learnt <- genes_tobe_learnt[genes_idx]
  mappedTSS <- mappedTSS[mappedTSS_idx,]

  # determine the dominant TSS
  cluster_domTSS_tab <- findDomTSS(mappedTSS = mappedTSS,
                                   clusters = clusters,
                                   ncores = ncore)
  cluster_domTSS_tab_gr <- makeGRangesFromDataFrame(df = cluster_domTSS_tab,
                                                    keep.extra.columns = TRUE)


  ##############
  ## learning step
  ################

  learnt_dist <- mclapply(1:length(genes_tobe_learnt), function(iii){

    tryCatch({

      temp_gene_tobe_learnt <- genes_tobe_learnt[iii]

      # grab gene info
      temp_strand <- as.character(strand(temp_gene_tobe_learnt))
      temp_annotated_TSS <- ifelse(temp_strand == "+",
                                   start(temp_gene_tobe_learnt),
                                   end(temp_gene_tobe_learnt))


      temp_clusters <- query_based_on_single_regionInput(single_region = temp_gene_tobe_learnt,
                                                         query_regions = cluster_domTSS_tab_gr)
      temp_clusters <- data.frame(temp_clusters)

      if(temp_strand == "+"){
        all_dist_to_be_learnt <- temp_clusters$tss - temp_annotated_TSS
        filtered1_dist_to_be_learnt <- all_dist_to_be_learnt[which(all_dist_to_be_learnt>learning_min_dist &
                                                                     all_dist_to_be_learnt < learning_max_dist)]
        filtered2_dist_to_be_learnt <- min(filtered1_dist_to_be_learnt)

        idx_tolearn <- which(all_dist_to_be_learnt == filtered2_dist_to_be_learnt)
        return(data.frame(learnt_dist = filtered2_dist_to_be_learnt,
                          counts_cluster = temp_clusters$counts_cluster[idx_tolearn],
                          width = temp_clusters$width[idx_tolearn],
                          counts = temp_clusters$counts[idx_tolearn]
        ))
      }else{
        all_dist_to_be_learnt <- temp_annotated_TSS - temp_clusters$tss
        filtered1_dist_to_be_learnt <- all_dist_to_be_learnt[which(all_dist_to_be_learnt>learning_min_dist &
                                                               all_dist_to_be_learnt < learning_max_dist)]
        suppressWarnings({ filtered2_dist_to_be_learnt <- min(filtered1_dist_to_be_learnt) })

        idx_tolearn <- which(all_dist_to_be_learnt == filtered2_dist_to_be_learnt)
        return(data.frame(learnt_dist = filtered2_dist_to_be_learnt,
                          counts_cluster = temp_clusters$counts_cluster[idx_tolearn],
                          width = temp_clusters$width[idx_tolearn],
                          counts = temp_clusters$counts[idx_tolearn]
                          ))
      }



    },
    error=function(e){

      return( data.frame(learnt_dist = Inf,
                         counts_cluster = Inf,
                         width = Inf,
                         counts = Inf
      ))

      })

  },
  mc.cores=ncore)

  learnt_dist <- do.call(rbind,learnt_dist)
  learnt_dist <- learnt_dist[which(is.finite(learnt_dist$learnt_dist)),]

  weights_to_address <- learnt_dist$counts_cluster
  learnt_dist_toadd <- learnt_dist$learnt_dist

  if(learning_methods=="median"){
    learnt_dist_out <- round(stats::median(learnt_dist_toadd))
  }else{
    learnt_dist_out <- round(sum(learnt_dist_toadd*weights_to_address/sum(weights_to_address)))
  }

  return(learnt_dist_out)

}
