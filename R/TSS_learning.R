
## to learn the TSS candidates and genuine TSS distance for single-end samples
TSS_learning <- function(region,
                         tx,
                         temp_result,
                         bam_path,
                         gene_extension,
                         read_assignmnet_dist,
                         read_number){
  chr <- seqnames(region)
  strand <- as.character(strand(region))
  start <- start(region)
  end <- end(region)

  # filtering peaks based on read number
  read_in <- reading_reads_loc(region = region,
                               bam_path = bam_path,
                               if_paired = FALSE,
                               gene_extension = gene_extension)
  read_in <- read_in[which(read_in$strand!=strand),]

  if (strand == "+"){
    temp_result$cellCounts <- lapply(temp_result$tss_candi,function(x){
      reads_distance <- abs(x - read_in$start)
      idx <- which(reads_distance <= read_assignmnet_dist)

      return(read_in$CB[idx])
    })
  }else{
    temp_result$cellCounts <- lapply(temp_result$tss_candi,function(x){
      reads_distance <- abs(x - read_in$end)
      idx <- which(reads_distance <= read_assignmnet_dist)
      return(read_in$CB[idx])
    })
  }
  read_in<- NULL
  temp_result$realReads <- lengths(temp_result$cellCounts)
  temp_result$realReads[which(is.na(temp_result$realReads))] <- 0
  temp_result$cellCounts <- NULL

  temp_result <- temp_result[which(temp_result$realReads>read_number),]

  # record reads distance
  if(nrow(temp_result)>0){
    if(strand=="+"){
      tt <- temp_result$tss_candi-unique(unlist(start(tx)))
    }else{
      tt <- unique(unlist(end(tx))) - temp_result$tss_candi
    }
    tt <- tt[which.min(abs(tt))]
    return(dist=tt)
  }
}
