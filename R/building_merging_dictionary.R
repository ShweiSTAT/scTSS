## build up a tss merging dictionary
building_merging_dictionary <- function(tss_names,tss_filtering_distance){
  tss_dictionary <- data.frame(do.call(rbind, strsplit(tss_names, ":", fixed=TRUE)))
  colnames(tss_dictionary) <- c("chr","gene_id","tss_candi_before","strand")
  tss_dictionary$tss_candi_before <- as.numeric(tss_dictionary$tss_candi_before)
  strand <- tss_dictionary$strand[1]
  if (strand == "-"){
    tss_dictionary <- tss_dictionary[order(tss_dictionary$tss_candi_before,decreasing = T),]
  }else{
    tss_dictionary <- tss_dictionary[order(tss_dictionary$tss_candi_before),]
  }

  ## make a interval bench
  int <- abs (tss_dictionary$tss_candi_before[-1] - tss_dictionary$tss_candi_before[-nrow(tss_dictionary)])

  # merge from the 5'- preparing merging dictionary
  if(length(int)==0){
    tss_dictionary$group <- 1
  }else{
    timer <- 1
    group.bentch <- 1
    dist.bentch <- 0
    group <- rep(0,nrow(tss_dictionary))
    group[1] <- group.bentch
    repeat{
      dist.bentch <- dist.bentch + int[timer]
      if (dist.bentch <=tss_filtering_distance){
        group[timer+1] <- group.bentch
        timer <- timer +1
      }else{
        group.bentch <- group.bentch +1
        group[timer+1] <- group.bentch
        timer <- timer +1
        dist.bentch <- 0
      }
      if(timer>= nrow(tss_dictionary)){break}
    }
    tss_dictionary$group <- group
  }

  tss_dictionary$tss_candi_start <- 0
  tss_dictionary$tss_candi_end <- 0
  group <- unique(tss_dictionary$group)
  for(y in group){
    temp <- tss_dictionary[which(tss_dictionary$group==y),]$tss_candi_before
    if (length(temp) >1){
      tss_dictionary[which(tss_dictionary$group==y),]$tss_candi_start <- min(temp)
      tss_dictionary[which(tss_dictionary$group==y),]$tss_candi_end <- max(temp)
    }
    if(length(temp) ==1){
      tss_dictionary[which(tss_dictionary$group==y),]$tss_candi_start <- temp
      tss_dictionary[which(tss_dictionary$group==y),]$tss_candi_end <- temp+1
    }
  }
  tss_dictionary <- data.frame(unmerged_tss= paste0(tss_dictionary[,1],":",tss_dictionary[,2],":",
                                                    tss_dictionary[,3],":",tss_dictionary[,4]),
                               merged_tss= paste0(tss_dictionary[,1],":",tss_dictionary[,2],":",
                                                  tss_dictionary[,6],":",tss_dictionary[,7],":",
                                                  tss_dictionary[,4]))
  return(tss_dictionary)
}
