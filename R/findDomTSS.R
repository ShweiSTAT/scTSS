# a function to find dominant TSS for clusters of a sample

findDomTSS <- function(mappedTSS,
                       clusters,
                       ncores = 1){
  # make GRanges
  gr <- makeGRangesFromDataFrame(df = clusters)
  tss_gr <- GRanges(seqnames = mappedTSS$seqnames,
                    ranges = IRanges(start = mappedTSS$tss, end = mappedTSS$tss),
                    strand = mappedTSS$strand)

  # Find overlaps between your GRanges object and the TSS GRanges
  overlaps <- findOverlaps(gr, tss_gr)
  overlaps <- data.frame(overlaps)
  cluster_idx <- unique(overlaps$queryHits)

  dom_tss <- mclapply(cluster_idx, function(x){
    tss_idx <- overlaps$subjectHits[which(overlaps$queryHits==x)]
    temp_INcluster_tss <- mappedTSS[tss_idx,]
    return(cbind(temp_INcluster_tss[which.max(temp_INcluster_tss[,4]),],
                 counts_cluster = sum(temp_INcluster_tss[,4])))
  },
  mc.cores = ncores
  )

  dom_tss <- do.call(rbind,dom_tss)
  dom_tss <- dom_tss[c(2,4,5)]
  dom_clusters <- clusters[cluster_idx,]
  outs <- cbind(dom_clusters,dom_tss)

  return(outs)
}
