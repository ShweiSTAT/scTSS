## to adjust the TSS based on the learnt distance for single-end samples
tss_adj_function <- function(distance_adjust,
                             gene_extension,
                             strand,
                             xn,
                             tss_candidates){
  # get the exon location
  exon_list <- data.frame(xn)
  if(strand=="+"){
    exon_list[1,]$start <- exon_list[1,]$start-gene_extension
  }else{
    exon_list[nrow(exon_list),]$end <- exon_list[nrow(exon_list),]$end+gene_extension
  }
  exon_list$width <- exon_list$end-exon_list$start
  # get tss candidate location
  for(ii in 1:length(tss_candidates)){
    # get the current location for tss candidates
    x <- tss_candidates[ii]
    if(strand=="+"){
      # initiate the distance to be adjusted
      dist_to_adjust <- distance_adjust
      exons_start_loc <- exon_list$start
      repeat{
        to_exon_start_dist <- x-exons_start_loc
        to_exon_start_dist <- to_exon_start_dist[which(to_exon_start_dist>0)]
        if(length(to_exon_start_dist)==0){
          tss_candidates[ii] <- exons_start_loc[1]
          break
        }
        min_exon_dist <- min(to_exon_start_dist)
        if(min_exon_dist>=dist_to_adjust){
          tss_candidates[ii] <- x-dist_to_adjust
          break
        }else{
          x <- exon_list$end[max(1,which(x-exons_start_loc==min_exon_dist)-1)]
          dist_to_adjust <- dist_to_adjust-min_exon_dist
        }
      }

    }else{
      # initiate the distance to be adjusted
      dist_to_adjust <- distance_adjust
      exons_start_loc <- exon_list$end
      repeat{
        to_exon_start_dist <- exons_start_loc - x
        to_exon_start_dist <- to_exon_start_dist[which(to_exon_start_dist>0)]
        if(length(to_exon_start_dist)==0){
          tss_candidates[ii] <- exons_start_loc[1]
          break
        }
        min_exon_dist <- min(to_exon_start_dist)
        if(min_exon_dist>=dist_to_adjust){
          tss_candidates[ii] <- x+dist_to_adjust
          break
        }else{
          x <- exon_list$start[min(length(exons_start_loc),which(exons_start_loc-x==min_exon_dist)+1)]
          dist_to_adjust <- dist_to_adjust-min_exon_dist
        }
      }

    }

  }
  return(tss_candidates)
}
