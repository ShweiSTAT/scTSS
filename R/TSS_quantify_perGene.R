# a function for quantify TSS clusters per Gene

TSS_quantify_perGene <- function(region,
                                 temp_result,
                                 reads_percentage,
                                 read_number,
                                 gene_extension,
                                 bam_path,
                                 cell_meta,
                                 if_paired){

  barcodes <- cell_meta$barcode
  if(if_paired){
    ## for paired-end data
    strand <- as.character(strand(region))
    temp_result <- data.frame(temp_result)
    read_in <- suppressWarnings({reading_reads_loc(region = region,
                                                   bam_path = bam_path,
                                                   if_paired = if_paired,
                                                   gene_extension = gene_extension)})

    if (strand == "+"){
      temp_result$cellCounts <- lapply(1:nrow(temp_result),function(ii){
        x <- temp_result[ii,]
        reads_distance <- sign(as.numeric(x[2]) - read_in$start.first)*sign(as.numeric(x[3]) - read_in$start.first)
        idx <- which(reads_distance < 0)
        return(read_in$CB.last[idx])
      })
    }else{
      temp_result$cellCounts <- lapply(1:nrow(temp_result),function(ii){
        x <- temp_result[ii,]
        reads_distance <- sign(as.numeric(x[2]) - read_in$end.last)*sign(as.numeric(x[3]) - read_in$end.last)
        idx <- which(reads_distance < 0)
        return(read_in$CB.first[idx])
      })
    }

    read_in<- NULL
    temp_result$realReads <- lengths(temp_result$cellCounts)
    temp_result$realReads[which(is.na(temp_result$realReads))] <- 0
    temp_result$PeakID <-  paste0(seqnames(region),":",
                                  names(region),":",
                                  temp_result$start, ":",
                                  temp_result$end,":",
                                  strand)


    # filtering based on reads number
    temp_result <- temp_result[which(temp_result$realReads>0),]
    temp_result <- temp_result[which(temp_result$realReads>max(reads_percentage*sum(temp_result$realReads),read_number)),]

    ## output tss candidates in a matrix format
    if(nrow(temp_result) > 0){
      geneOuts <- data.frame(matrix(ncol = length(barcodes)+1, nrow = nrow(temp_result)))
      colnames(geneOuts) <- c("TSS_clusters",barcodes)
      geneOuts$TSS_clusters <- temp_result$PeakID
      for(ii in 1:nrow(geneOuts)){
        temp_table <- table(temp_result$cellCounts[ii])
        geneOuts[ii,2:(length(barcodes)+1)] <- temp_table[match(colnames(geneOuts)[-1],names(temp_table))]
        temp_table <- NULL
      }
      setDT(geneOuts,key = "TSS_clusters")
      temp_result <- NULL
      removeNAs(geneOuts)
      return(geneOuts)
    }else{
      stop()
    }
  }else{
    ## for single-end data
    strand <- as.character(strand(region))
    temp_result <- data.frame(temp_result)
    read_in <- suppressWarnings({reading_reads_loc(region = region,
                                                   bam_path = bam_path,
                                                   if_paired = FALSE,
                                                   gene_extension = gene_extension)})

    if (strand == "+"){
      temp_result$cellCounts <- lapply(1:nrow(temp_result),function(ii){
        x<- temp_result[ii,]
        reads_distance <- sign(as.numeric(x[2]) - read_in$start)*sign(as.numeric(x[3]) - read_in$start)
        idx <- which(reads_distance < 0)
        return(read_in$CB[idx])
      })
    }else{
      temp_result$cellCounts <- lapply(1:nrow(temp_result),function(ii){
        x<- temp_result[ii,]
        reads_distance <- sign(as.numeric(x[2]) - read_in$end)*sign(as.numeric(x[3]) - read_in$end)
        idx <- which(reads_distance < 0)
        return(read_in$CB[idx])
      })
    }

    read_in<- NULL
    temp_result$realReads <- lengths(temp_result$cellCounts)
    temp_result$realReads[which(is.na(temp_result$realReads))] <- 0
    temp_result$PeakID <-  paste0(as.character(seqnames(region)),":",
                                  names(region),":",
                                  temp_result$tss5prime, ":",
                                  temp_result$width,":",
                                  strand)

    temp_result <- temp_result[which(temp_result$realReads>0),]
    # # filtering based on reads number
    temp_result <- temp_result[which(temp_result$realReads>max(reads_percentage*sum(temp_result$realReads),read_number)),]


    ## output tss candidates in a matrix format
    if(nrow(temp_result) > 0){
      geneOuts <- data.frame(matrix(ncol = length(barcodes)+1, nrow = nrow(temp_result)))
      colnames(geneOuts) <- c("TSS_clusters",barcodes)
      geneOuts$TSS_clusters <- temp_result$PeakID
      for(ii in 1:nrow(geneOuts)){
        temp_table <- table(temp_result$cellCounts[ii])
        geneOuts[ii,2:(length(barcodes)+1)] <- temp_table[match(colnames(geneOuts)[-1],names(temp_table))]
        temp_table <- NULL
      }
      setDT(geneOuts,key = "TSS_clusters")
      removeNAs(geneOuts)
      temp_result <- NULL

      return(geneOuts)
    }else{
      stop()
    }

  }

}
