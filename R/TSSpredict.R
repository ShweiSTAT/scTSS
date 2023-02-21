#' Splicing-aware TSS prediction from .bam file
#'
#' @param read_assignmnet_dist An integer specifying the threshold distance we
#' use to assign reads to TSS candidates. Defaults to 200.
#' @param gene_extension An integer specifying the extension before the
#' annotated 5' of a gene. Defaults to 500.
#' @param tss_filtering_distance An integer specifying the threshold distance
#' for TSS candidates that should be considered from the same genuine TSS.
#' This should be set to the read length of the sequence technology.
#' Defaults to 150.
#' @param reads_percentage A percentage specifying the minim percentage of
#' reads on a gene should be assigned to a qualified TSS candidates.
#' Defaults to 0.05.
#' @param read_number An integer specifying the minim number of reads that
#' should be assigned to a qualified TSS candidate. Defaults to 50.
#' @param gene_of_interest A character or a vector of characters, specifying
#' the genes that are of interest in terms of predicting TSS.
#' The default is NULL, which means we need to check all genes provided by the
#' the genome annotation (see anno_path). Otherwise provide a subset of gene IDs
#' that is in the provided annotation.
#' @param if_paired FALSE/TRUE, specifying if we have paired-end data (TRUE)
#' or single-end data (FALSE). Default set to TRUE.
#' @param ncore An integer specifying the number of cores used. Defaults to 1.
#' @param bam_path The path to the .bam, which should be stored under the same
#' directory as its .bai file.
#' @param tss_path The path to the .txt file that stores the information
#' of TSS candidates
#' @param anno_path The path to the genome annotation file, which should be in
#' the .gtf format.
#' @param barcodes_path The path to the barcode file, which should be a .csv
#' file with its first column being the barcode and no header.
#' @return \code{TSSprediction} outputs a matrix of read counts for the
#' predicted TSS by cell. This matrix is stored in data.table format.
#' @export TSSpredict
#' @importFrom data.table  .SD := data.table rbindlist set setDT
#' @import parallel
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom Rsamtools ScanBamParam
#' @importFrom utils read.csv read.delim head tail
#' @importFrom stats median
#'
#' @author Shiwei Fu, \email{shiwei.fu@email.ucr.edu}



TSSpredict <- function(read_assignmnet_dist = 200,
                           gene_extension = 500,
                           tss_filtering_distance = 150,
                           reads_percentage = 0.05,
                           read_number = 50,
                           gene_of_interest = NULL,
                           if_paired = TRUE,
                           ncore = 1,
                           bam_path,
                           anno_path,
                           tss_path,
                           barcodes_path){
  ## load in the needed files from path
  txdb <- makeTxDbFromGFF(anno_path)
  # get gene info
  gene_list <- genes(txdb)
  if(is.null(gene_of_interest)){
    all_genes <- names(gene_list)
  }else{
    if(sum(!gene_of_interest%in%names(gene_list))==0){
      all_genes <- gene_of_interest
    }else{
      stop("ERROR: gene_of_interest is not valide. Please provide real gene names.")
    }
    }
  # extract exon information
  exons_list_by_gene <- reduce(exonsBy(txdb, by="gene"))
  # extract_tx information
  tx_list_by_gene <- transcriptsBy(txdb,by=c("gene"))
  # outs from Homer
  results <- read.delim(tss_path,sep="\t")
  # loading barcodes
  barcodes <- read.csv(barcodes_path,header = F)[,1]
  # loading one-TSS gene
  if(!if_paired){
    print("This is a single-end data. We are extracting the one-TSS genes for TSS learning...")
    oneTSSgene <- getOneTSSgenes(tx_list_by_gene = tx_list_by_gene,
                                 ncore = ncore)
  }

  if(if_paired){
    print("This is a paired-end data. We are predicting TSS and counting reads...")
  tss_outs <- mclapply(1:length(all_genes),function(iii){
    tryCatch({
      region <- gene_list[all_genes[iii]]
      chr <- as.character(seqnames(region))
      strand <- as.character(strand(region))
      start <- start(region)
      end <- end(region)

      tx <- tx_list_by_gene[all_genes[iii]]

      xn <- exons_list_by_gene[all_genes[iii]]

      temp_result <- results[which(results$chr==chr),]
      if(strand == "+"){
        temp_result <- temp_result[which((temp_result$start>start-gene_extension) & (temp_result$end < end)),]
      }else{
        temp_result <- temp_result[which((temp_result$start>start) & (temp_result$end < end+gene_extension)),]
      }
      temp_result <- temp_result[order(temp_result$start),]
      temp_result$tss_candi <- round((temp_result$start+temp_result$end)/2)

      outs <- tss_candidate_handle(region = region,
                                   tx = tx,
                                   xn = xn,
                                   temp_result = temp_result,
                                   read_assignmnet_dist = read_assignmnet_dist,
                                   tss_filtering_distance = tss_filtering_distance,
                                   gene_extension = gene_extension,
                                   reads_percentage = reads_percentage,
                                   read_number = read_number,
                                   distance_adjust = NULL,
                                   if_paired = if_paired,
                                   bam_path = bam_path,
                                   barcodes = barcodes)
      if (iii%%500 == 0){
        gc()
      }
      return(outs)
    },
    error=function(e){cat("Gene: ",all_genes[iii]," excluded", "\n")})
  },mc.cores = ncore)
  tss_outs <- tss_outs[lengths(tss_outs) >0]
  tss_outs <- rbindlist(tss_outs,use.names=T)
  setDT(tss_outs,"unmerged_tss")

  }else{
    ## 1) The learning step
    print("We are learning the adjustment distance for your input single-end data...")
    tss_learnt <- mclapply(1:length(oneTSSgene),function(iii){
      tryCatch({
        region <- gene_list[oneTSSgene[iii]]
        chr <- as.character(seqnames(region))
        strand <- as.character(strand(region))
        start <- start(region)
        end <- end(region)

        tx <- tx_list_by_gene[oneTSSgene[iii]]

        xn <- exons_list_by_gene[oneTSSgene[iii]]

        temp_result <- results[which(results$chr==chr),]
        if(strand == "+"){
          temp_result <- temp_result[which((temp_result$start>start-gene_extension) & (temp_result$end < end)),]
        }else{
          temp_result <- temp_result[which((temp_result$start>start) & (temp_result$end < end+gene_extension)),]
        }
        temp_result <- temp_result[order(temp_result$start),]
        temp_result$tss_candi <- round((temp_result$start+temp_result$end)/2)

        outs <- TSS_learning(region = region,
                             tx = tx,
                             temp_result = temp_result,
                             bam_path = bam_path,
                             gene_extension = gene_extension,
                             read_assignmnet_dist = read_assignmnet_dist,
                             read_number = read_number)
        if (iii%%500 == 0){
          gc()
        }
        return(outs)
      },
      error=function(e){cat("Gene: ",oneTSSgene[iii]," is not used in the learning step", "\n")})
    },mc.cores = ncore)
    tss_learnt <- tss_learnt[lengths(tss_learnt)>0]
    tss_learnt <- unlist(tss_learnt)
    tss_learnt <- tss_learnt[which(tss_learnt>0 & abs(tss_learnt)<600)]
    tss_learnt <- ceiling(median(tss_learnt))
    print(paste0("The learnt adjustment distance for sample ",
                 bam_path," is ",
                 tss_learnt, "bp."))


    ## 2) The predicting step
    print("We are predicting TSS and counting reads...")
    tss_outs <- mclapply(1:length(all_genes),function(iii){
      tryCatch({
        region <- gene_list[all_genes[iii]]
        chr <- as.character(seqnames(region))
        strand <- as.character(strand(region))
        start <- start(region)
        end <- end(region)

        tx <- tx_list_by_gene[all_genes[iii]]

        xn <- exons_list_by_gene[all_genes[iii]]

        temp_result <- results[which(results$chr==chr),]
        if(strand == "+"){
          temp_result <- temp_result[which((temp_result$start>start-gene_extension) & (temp_result$end < end)),]
        }else{
          temp_result <- temp_result[which((temp_result$start>start) & (temp_result$end < end+gene_extension)),]
        }
        temp_result <- temp_result[order(temp_result$start),]
        temp_result$tss_candi <- round((temp_result$start+temp_result$end)/2)

        outs <- tss_candidate_handle(region = region,
                                     tx = tx,
                                     xn = xn,
                                     temp_result = temp_result,
                                     read_assignmnet_dist = read_assignmnet_dist,
                                     tss_filtering_distance = tss_filtering_distance,
                                     gene_extension = gene_extension,
                                     reads_percentage = reads_percentage,
                                     read_number = read_number,
                                     distance_adjust = tss_learnt,
                                     if_paired = if_paired,
                                     bam_path = bam_path,
                                     barcodes = barcodes)
        if (iii%%500 == 0){
          gc()
        }
        return(outs)
      },
      error=function(e){cat("Gene: ",all_genes[iii]," excluded", "\n")})
    },mc.cores = ncore)
    tss_outs <- tss_outs[lengths(tss_outs) >0]
    tss_outs <- rbindlist(tss_outs,use.names=T)
    setDT(tss_outs,"unmerged_tss")
  }

  print("TSS prediction finished!")
  return(tss_outs)
}

