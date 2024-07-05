# a function to adjust near-site TSS clusters for a single gene

TSS_adjust_perGene <- function(this_gene_region,
                               this_gene_clusters,
                               this_gene_xn,
                               this_sample_leart_dist,
                               gene_extension = 500){

  strand <- as.character(strand(this_gene_region))
  temp_result <- as.data.frame(this_gene_clusters)

  # extend the first exon, in case there are intergenetic TSSs
  exon_list <- data.frame(this_gene_xn)
  if(strand=="+"){
    exon_list[1,]$start <- exon_list[1,]$start-gene_extension
  }else{
    exon_list[nrow(exon_list),]$end <- exon_list[nrow(exon_list),]$end+gene_extension
  }
  exon_list$width <- exon_list$end-exon_list$start+1


  #### step 1: decide if the start site of the cluster is on the exon; IF not remove the cluster

  for(ii in 1:nrow(temp_result)){
    temp <- temp_result[ii,]

    if(strand == "+"){
      start_xn_info <- check_loc_exons(exons = exon_list,
                                       loc = temp$start)
    }else{
      start_xn_info <- check_loc_exons(exons = exon_list,
                                       loc = temp$end)
    }

    if(!start_xn_info$in_exon){
      temp_result <- temp_result[-ii,]
    }

  }


  #### step 2: move the start site of the genomic regions along the exons
  if(strand == "+"){
    start_sites <- temp_result$start
  }else{
    start_sites <- temp_result$end
  }
  org_clusters_location <- temp_result


  for(ii in 1:length(start_sites)){
    # get the current location for tss candidates
    x <- start_sites[ii]
    if(strand=="+"){
      # initiate the distance to be adjusted
      dist_to_adjust <- this_sample_leart_dist
      exons_start_loc <- exon_list$start
      repeat{
        to_exon_start_dist <- x-exons_start_loc
        to_exon_start_dist <- to_exon_start_dist[which(to_exon_start_dist>0)]
        if(length(to_exon_start_dist)==0){
          start_sites[ii] <- exons_start_loc[1]
          break
        }
        min_exon_dist <- min(to_exon_start_dist)
        if(min_exon_dist>=dist_to_adjust){
          start_sites[ii] <- x-dist_to_adjust
          break
        }else{
          x <- exon_list$end[max(1,which(x-exons_start_loc==min_exon_dist)-1)]
          dist_to_adjust <- dist_to_adjust-min_exon_dist
        }
      }

    }else{
      # initiate the distance to be adjusted
      dist_to_adjust <- this_sample_leart_dist
      exons_start_loc <- exon_list$end
      repeat{
        to_exon_start_dist <- exons_start_loc - x
        to_exon_start_dist <- to_exon_start_dist[which(to_exon_start_dist>0)]
        if(length(to_exon_start_dist)==0){
          start_sites[ii] <- exons_start_loc[length(exons_start_loc)]
          break
        }
        min_exon_dist <- min(to_exon_start_dist)
        if(min_exon_dist>=dist_to_adjust){
          start_sites[ii] <- x+dist_to_adjust
          break
        }else{
          x <- exon_list$start[min(length(exons_start_loc),which(exons_start_loc-x==min_exon_dist)+1)]
          dist_to_adjust <- dist_to_adjust-min_exon_dist
        }
      }

    }

  }

  #### step 3: reconstruct clusters
  if(strand == "+"){
    org_start_sites <- org_clusters_location$start
  }else{
    org_start_sites <- org_clusters_location$end
  }

  dist_change <- org_start_sites - start_sites

  out <- data.frame(seqnames = temp_result$seqnames[1],
                    tss5prime = start_sites,
                    width = temp_result$width,
                    strand = strand,
                    org_start  = temp_result$start,
                    org_end = temp_result$end)

  return(out)
}
