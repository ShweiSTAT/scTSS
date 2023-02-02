
# a peak selection function
tss_candidate_handle <- function(region,
                                 tx,
                                 xn,
                                 temp_result,
                               read_assignmnet_dist,
                               tss_filtering_distance,
                               gene_extension,
                               reads_percentage,
                               read_number,
                               distance_adjust = NULL,
                               if_paired,
                               bam_path,
                               barcodes){
  if(if_paired){
    chr <- as.character(seqnames(region))
    strand <- as.character(strand(region))
    start <- start(region)
    end <- end(region)

    ## prestep: assign reads to all tss candidates; renames tss id
    read_in <-  reading_reads_loc(region = region,
                                  bam_path = bam_path,
                                  if_paired = if_paired,
                                  gene_extension = gene_extension)
    read_in <- read_in[which(read_in$strand.first==strand),]
    if (strand == "+"){
      temp_result$cellCounts <- lapply(temp_result$tss_candi,function(x){
        reads_distance <- abs(x - read_in$start.first)
        idx <- which(reads_distance <= read_assignmnet_dist)

        return(read_in$CB.last[idx])
      })
    }else{
      temp_result$cellCounts <- lapply(temp_result$tss_candi,function(x){
        reads_distance <- abs(x - read_in$end.last)
        idx <- which(reads_distance <= read_assignmnet_dist)
        return(read_in$CB.last[idx])
      })
    }
    read_in<- NULL
    temp_result$realReads <- lengths(temp_result$cellCounts)
    temp_result$realReads[which(is.na(temp_result$realReads))] <- 0
    temp_result$PeakID <-  paste0(chr,":",names(region),":",temp_result$tss_candi,":",strand)

    temp_result <- temp_result[which(temp_result$realReads>0),]

    ## peak selection
    # step1: extend the first exon
    exon_list <- data.frame(xn)
    if(strand=="+"){
      exon_list[1,]$start <- exon_list[1,]$start-gene_extension
    }else{
      exon_list[nrow(exon_list),]$end <- exon_list[nrow(exon_list),]$end+gene_extension
    }
    exon_list$width <- exon_list$end-exon_list$start


    # step2: merge within-range intronic candidates to the nearest exon edge and mark them exonic
    tss.candidates <- temp_result$tss_candi

    exon_marker <- rep(FALSE,nrow(temp_result))
    names(exon_marker) <- temp_result$PeakID
    ex <- c(exon_list$start,exon_list$end)

    for(ii in 1:nrow(temp_result)){
      x <- tss.candidates[ii]
      if (sum(sign(exon_list$start - x)*sign(exon_list$end - x)<0) >0){
        exon_marker[ii] <- TRUE
      }else{
        distance_to_exon <- abs(ex - x)
        replace.idx <- which.min(distance_to_exon)
        if(distance_to_exon[replace.idx]<tss_filtering_distance){
          tss.candidates[ii] <- ex[replace.idx]
          exon_marker[ii] <- TRUE
        }
      }
    }
    temp_result$exon_marker<-exon_marker
    temp_result$tss_candi <- tss.candidates


    # step3: peak filtration on tx level(on-exon tss candidates)/ genome level(off exon tss candidates)
    exon_filter <- as.numeric()
    intron_filter <- as.numeric()

    # --step 3a: get back to tx loc for on-exon candidates
    for (ii in 1: nrow(temp_result)){
      temp_peak <- temp_result[ii,]
      tss.candidates <- temp_peak$tss_candi
      if_on_exon <- temp_peak$exon_marker

      if (strand == "+" & if_on_exon){
        exon_idx <- which(exon_list$start - tss.candidates <= 0)
        tx.loc <- sum(exon_list$width[head(exon_idx,-1)]) + tss.candidates - exon_list$start[tail(exon_idx,1)]
        exon_filter <- append(exon_filter,tx.loc)
      }
      if(strand == "-" & if_on_exon){
        exon_idx <- which(exon_list$end - tss.candidates >= 0)
        tx.loc <- sum(exon_list$width[tail(exon_idx,-1)]) + exon_list$end[head(exon_idx,1)] - tss.candidates
        exon_filter <- append(exon_filter,tx.loc)
      }
    }

    exon_filter <- data.frame(filter=exon_filter,
                              reads = temp_result$realReads[which(exon_marker==TRUE)],
                              row.names =temp_result$PeakID[which(exon_marker==TRUE)])
    exon_filter$cellCounts <- temp_result[which(exon_marker==TRUE),]$cellCounts
    exon_filter <- exon_filter[order(exon_filter$filter),]

    intron_filter <- temp_result$tss_candi[which(exon_marker==FALSE)]
    intron_filter <- data.frame(filter=intron_filter,
                                reads=temp_result$realReads[which(exon_marker==FALSE)],
                                row.names = temp_result$PeakID[which(exon_marker==FALSE)])
    intron_filter$cellCounts <- temp_result[which(exon_marker==FALSE),]$cellCounts
    if(strand=="+"){
      intron_filter <- intron_filter[order(intron_filter$filter),]
    }else{
      intron_filter <- intron_filter[order(intron_filter$filter,decreasing = T),]
    }
    temp_result$cellCounts <- NULL

    # --step 3b: filtering peaks
    exon_filter <- peak_filter(exon_filter,tss_filtering_distance)
    intron_filter <- peak_filter(intron_filter,tss_filtering_distance)
    filter <- rbind(exon_filter,intron_filter)

    tss_loss_selection <- nrow(temp_result) - length(c(rownames(exon_filter),rownames(intron_filter)))
    exon_filter <- NULL
    intron_filter <- NULL

    temp_result <- temp_result[which(temp_result$PeakID%in%rownames(filter)),]
    temp_result$cellCounts <- filter[match(temp_result$PeakID,rownames(filter)),]$cellCounts
    temp_result <- temp_result[order(temp_result$start),]

    # step4: filtering based on reads number
    tss_loss_filtering <- nrow(temp_result)-length(which(temp_result$realReads>max(reads_percentage*sum(temp_result$realReads),read_number)))
    temp_result <- temp_result[which(temp_result$realReads>max(reads_percentage*sum(temp_result$realReads),read_number)),]

    ## output tss candidates in a matrix format
    geneOuts <- data.frame(matrix(ncol = length(barcodes)+1, nrow = nrow(temp_result)))
    colnames(geneOuts) <- c("unmerged_tss",barcodes)
    geneOuts$unmerged_tss <- temp_result$PeakID
    for(ii in 1:nrow(geneOuts)){
      temp_table <- table(temp_result[ii,]$cellCounts)
      geneOuts[ii,2:(length(barcodes)+1)] <- temp_table[match(colnames(geneOuts)[-1],names(temp_table))]
      temp_table <- NULL
    }
    setDT(geneOuts,key = "unmerged_tss")
    removeNAs(geneOuts)
    temp_result <- NULL

    return(geneOuts)
  }else{
    chr <- as.character(seqnames(region))
    strand <- as.character(strand(region))
    start <- start(region)
    end <- end(region)

    ## prestep: assign reads to all TSS candidates; renames TSS id
    read_in <- reading_reads_loc(region = region,
                                 bam_path = bam_path,
                                 if_paired = if_paired,
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
    temp_result$PeakID <-  paste0(chr,":",names(region),":",temp_result$tss_candi,":",strand)

    temp_result <- temp_result[which(temp_result$realReads>0),]

    # splicing-aware and false positive controlled tss selection
    # step1: extend the first exon
    exon_list <- data.frame(xn)
    if(strand=="+"){
      exon_list[1,]$start <- exon_list[1,]$start-gene_extension
    }else{
      exon_list[nrow(exon_list),]$end <- exon_list[nrow(exon_list),]$end+gene_extension
    }
    exon_list$width <- exon_list$end-exon_list$start


    # step2: merge within-range intronic candidates to the nearest exon edge and mark them exonic
    tss.candidates <- temp_result$tss_candi

    exon_marker <- rep(FALSE,nrow(temp_result))
    names(exon_marker) <- temp_result$PeakID
    ex <- c(exon_list$start,exon_list$end)

    for(ii in 1:nrow(temp_result)){
      x <- tss.candidates[ii]
      if (sum(sign(exon_list$start - x)*sign(exon_list$end - x)<0) >0){
        exon_marker[ii] <- TRUE
      }else{
        distance_to_exon <- abs(ex - x)
        replace.idx <- which.min(distance_to_exon)
        if(distance_to_exon[replace.idx]<tss_filtering_distance){
          tss.candidates[ii] <- ex[replace.idx]
          exon_marker[ii] <- TRUE
        }
      }
    }
    temp_result$exon_marker<-exon_marker
    temp_result$tss_candi <- tss.candidates


    # step3: peak filtration on tx level(on-exon tss candidates)/ genome level(off exon tss candidates)
    exon_filter <- as.numeric()
    intron_filter <- as.numeric()

    # --step 3a: get back to tx loc for on-exon candidates
    for (ii in 1: nrow(temp_result)){
      temp_peak <- temp_result[ii,]
      tss.candidates <- temp_peak$tss_candi
      if_on_exon <- temp_peak$exon_marker

      if (strand == "+" & if_on_exon){
        exon_idx <- which(exon_list$start - tss.candidates <= 0)
        tx.loc <- sum(exon_list$width[head(exon_idx,-1)]) + tss.candidates - exon_list$start[tail(exon_idx,1)]
        exon_filter <- append(exon_filter,tx.loc)
      }
      if(strand == "-" & if_on_exon){
        exon_idx <- which(exon_list$end - tss.candidates >= 0)
        tx.loc <- sum(exon_list$width[tail(exon_idx,-1)]) + exon_list$end[head(exon_idx,1)] - tss.candidates
        exon_filter <- append(exon_filter,tx.loc)
      }
    }

    exon_filter <- data.frame(filter=exon_filter,
                              reads = temp_result$realReads[which(exon_marker==TRUE)],
                              row.names =temp_result$PeakID[which(exon_marker==TRUE)])
    exon_filter$cellCounts <- temp_result[which(exon_marker==TRUE),]$cellCounts
    exon_filter <- exon_filter[order(exon_filter$filter),]

    intron_filter <- temp_result$tss_candi[which(exon_marker==FALSE)]
    intron_filter <- data.frame(filter=intron_filter,
                                reads=temp_result$realReads[which(exon_marker==FALSE)],
                                row.names = temp_result$PeakID[which(exon_marker==FALSE)])
    intron_filter$cellCounts <- temp_result[which(exon_marker==FALSE),]$cellCounts
    if(strand=="+"){
      intron_filter <- intron_filter[order(intron_filter$filter),]
    }else{
      intron_filter <- intron_filter[order(intron_filter$filter,decreasing = T),]
    }
    temp_result$cellCounts <- NULL

    #--step 3b: filtering peaks
    exon_filter <- peak_filter(exon_filter,tss_filtering_distance)
    intron_filter <- peak_filter(intron_filter,tss_filtering_distance)
    filter <- rbind(exon_filter,intron_filter)

    tss_loss_selection <- nrow(temp_result) - length(c(rownames(exon_filter),rownames(intron_filter)))
    exon_filter <- NULL
    intron_filter <- NULL

    temp_result <- temp_result[which(temp_result$PeakID%in%rownames(filter)),]
    temp_result$cellCounts <- filter[match(temp_result$PeakID,rownames(filter)),]$cellCounts
    temp_result <- temp_result[order(temp_result$start),]

    # step4: filtering based on reads number
    tss_loss_filtering <- nrow(temp_result)-length(which(temp_result$realReads>max(reads_percentage*sum(temp_result$realReads),read_number)))
    temp_result <- temp_result[which(temp_result$realReads>max(reads_percentage*sum(temp_result$realReads),read_number)),]

    # final step: adjust TSS by the learnt distance
    temp_result$tss_candi <- tss_adj_function(distance_adjust = distance_adjust,
                                              gene_extension = gene_extension,
                                              strand = strand,
                                              xn=xn,
                                              tss_candidates = temp_result$tss_candi)
    temp_result$PeakID <-  paste0(chr,":",names(region),":",temp_result$tss_candi,":",strand)

    ## output TSS candidates in a matrix format
    geneOuts <- data.frame(matrix(ncol = length(barcodes)+1, nrow = nrow(temp_result)))
    colnames(geneOuts) <- c("unmerged_tss",barcodes)
    geneOuts$unmerged_tss <- temp_result$PeakID
    for(ii in 1:nrow(geneOuts)){
      temp_table <- table(temp_result[ii,]$cellCounts)
      geneOuts[ii,2:(length(barcodes)+1)] <- temp_table[match(colnames(geneOuts)[-1],names(temp_table))]
      temp_table <- NULL
    }
    setDT(geneOuts,key = "unmerged_tss")
    removeNAs(geneOuts)
    temp_result <- NULL

    return(geneOuts)
  }
}
