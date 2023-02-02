# to filter TSS based on distance and read counts
peak_filter <- function(filter,tss_filtering_distance){
  tx_loc <- filter[,"filter"]
  names(tx_loc) <- row.names(filter)
  reads <- filter[,"cellCounts"]
  timer <-1
  repeat{
    idx <- which(abs(tx_loc - tx_loc[timer])<= tss_filtering_distance)
    if (length(idx) >1){
      temp_reads <- reads[idx]
      reserve_idx <- idx[which.max(lengths(temp_reads))]
      reads[[reserve_idx]] <- Reduce(c,temp_reads)
      idx <- setdiff(idx,reserve_idx)
      tx_loc <- tx_loc[-idx]
      reads <- reads[-idx]
    }
    timer <- timer +1
    if(timer > length(tx_loc)){
      break
    }
  }
  outs <- data.frame(filter=tx_loc)
  outs$cellCounts <- reads
  return(outs)
}
