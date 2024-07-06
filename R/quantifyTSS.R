#' Quantify TSS clusters' expression levels for each cell
#'
#' @param anno_path The path to the genome annotation file, which should be in
#' the .gtf format.
#' @param tss_clusters A data.frame specifies the TSS clusters to be quantified
#' for this sample. For paired-end data, these columns should be included:
#' seqnames, start, end and strand. For single-end data, these column should be
#' included: seqnames, tss5prime, width, strand, org_start, org_end.
#' @param bam_path The path to the .bam, which should be stored under the same
#' directory as its .bai file.
#' @param cell_meta a data.frame specifies the meta information for each cell.
#' "barcode" and "sampleID" must be included as columns.
#' "barcode" specifies the barcodes for each of the cells being sequenced
#' in this sample; "sampleID" is a unique ID for each sample.
#' Other information are also suggested to be included
#' in this `data.frame` as well, such as cell types or the
#' sample's other information (age, gender, etc.).
#' @param gene_of_interest A character or a vector of characters, specifying
#' the gene names that are of interests. The input gene names should be same
#' as those in the genome annotation (see anno_path).
#' The default is NULL, which means we need to check all genes provided by the
#' genome annotation (see anno_path).
#' @param reads_percentage A percentage specifying the minim percentage of
#' reads on a gene that should be assigned to a valid TSS clusters within this sample.
#' This is used to filter out lowly expressed TSS clusters.
#' Defaults to 0.05.
#' @param read_number An integer specifying the minim number of reads that
#' should be assigned to a valid TSS clusters within this sample.
#' This is used to filter out lowly expressed TSS clusters. Defaults to 50.
#' @param gene_extension An integer specifying the length of the extended region
#' before the annotated 5' of a gene when quantifying TSS. Defaults to 500 (bp).
#' @param if_paired FALSE/TRUE, specifying if we have paired-end data (TRUE)
#' or single-end data (FALSE). Default set to TRUE.
#' @param ncore An integer specifying the number of cores used for computation.
#' Defaults to 1.
#'
#' @return \code{quantifyTSS} returns a list containing these two elements:
#' 1) TSS_cluster_counts, a data.table for the expression (counts) of each
#' TSS clusters at cell level. 2) Column_meta, a data.frame
#' containing the meta information for each of the columns (cells)
#' in TSS_cluster_counts. The columns in TSS_cluster_counts and
#' the rows in Column_meta have one on one correspondence.
#' @export quantifyTSS
#'
#' @import parallel
#' @import GenomicAlignments
#' @import GenomicRanges
#' @import GenomicFeatures
#' @importFrom data.table  .SD := data.table rbindlist set setDT setnames
#' @importFrom Rsamtools ScanBamParam
#'
quantifyTSS <- function(anno_path,
                        tss_clusters,
                        bam_path,
                        cell_meta,
                        gene_of_interest = NULL,
                        reads_percentage = 0.05,
                        read_number = 10,
                        gene_extension = 500,
                        if_paired = TRUE,
                        ncore = 1){


  if(if_paired){
    #### paired-end data
    ## load in data
    txdb <- makeTxDbFromGFF(anno_path)
    gene_list <- genes(txdb)

    # check if tss_clusters is loaded right
    if(!(sum(c("chr","start","end","strand")
             %in%colnames(tss_clusters)) == 4 |
         sum(c("seqnames","start","end","strand")
             %in%colnames(tss_clusters)) == 4 )){
      stop("\"tss_clusters\" has to be a data.frame with at least these columns: chr/seqnames, start, end, strand.")
    }


    # check if cell_meta is loaded right
    if(sum(!c("barcode","sampleID")%in%colnames(cell_meta))>1){
      stop("\"barcode\" and \"sampleID\" have to be included in arguement \"cell_meta\" as columns.")
    }


    # check if gene_of_interest is correct
    if(is.null(gene_of_interest)){
      all_genes <- names(gene_list)
    }else{
      if(sum(!gene_of_interest%in%names(gene_list))==0){
        all_genes <- gene_of_interest
        gene_list <- gene_list[all_genes]
      }else{
        stop(" \"gene_of_interest\" is not valide. Please provide real gene names.")
      }
    }

    ################### pre-step: reduce the genes according to the TSS clusters
    clusters_gr <- makeGRangesFromDataFrame(df = tss_clusters,
                                            keep.extra.columns = TRUE)
    resized_gene_list <- resize(gene_list,width(gene_list)+gene_extension,fix = "end")
    overlaps <- findOverlaps(resized_gene_list, clusters_gr)
    overlaps <- data.frame(overlaps)
    genes_idx <- unique(overlaps$queryHits)
    cluster_idx <- unique(overlaps$subjectHits)

    gene_list <- gene_list[genes_idx]
    tss_clusters <- tss_clusters[cluster_idx,]

    ## out reduced results: reduced genes
    all_genes <- names(gene_list)
    ## out reduced results：tss_clsuters
    clusters_gr <- makeGRangesFromDataFrame(df = tss_clusters,
                                            keep.extra.columns = TRUE)
    ######################## quantify: gene by gene
    print(paste0("Start quantifing ",length(clusters_gr)," clusters on ", length(gene_list)," genes"))
    tss_outs <- parallel:: mclapply(1:length(gene_list), function(iii){
      tryCatch({
      region <- gene_list[iii]
      temp_result <- query_based_on_single_regionInput(single_region = region,
                                                       query_regions = clusters_gr)
      outs <- TSS_quantify_perGene(region = region,
                                   temp_result = temp_result,
                                   reads_percentage = reads_percentage,
                                   read_number = read_number,
                                   gene_extension = gene_extension,
                                   bam_path = bam_path,
                                   cell_meta = cell_meta,
                                   if_paired = if_paired)
      if (iii%%500 == 0){
        gc()
      }

      return(outs)
      }, error=function(e){})

    },mc.cores = ncore)

    names(tss_outs) <- all_genes
    print("Counting reads finished!")
    tss_outs <- tss_outs[lengths(tss_outs)>0]
    tss_outs <- rbindlist(tss_outs,use.names=T)
    setDT(tss_outs,"TSS_clusters")

    setnames(tss_outs,
             cell_meta$barcode,
             paste0(cell_meta$sampleID,"_",cell_meta$barcode) )
    cell_meta$barcode <- paste0(cell_meta$sampleID,"_",cell_meta$barcode)

    return(list(TSS_clsuter_counts = tss_outs,
                Column_meta = cell_meta))

  }else{
    #### single-end data
    ## load in data
    txdb <- makeTxDbFromGFF(anno_path)
    gene_list <- genes(txdb)

    # check if tss_clusters is loaded right
    if(!(sum(c("seqnames","tss5prime","width","strand","org_start", "org_end")
             %in%colnames(tss_clusters)) == 6 |
         sum(c("chr","tss5prime","width","strand","org_start", "org_end")
             %in%colnames(tss_clusters)) == 6)){
      stop("\"tss_clusters\" has to be a data.frame with at least these columns: chr/seqnames, tss5prime, width, strand, org_start, org_end.")
    }


    # check if cell_meta is loaded right
    if(sum(!c("barcode","sampleID")%in%colnames(cell_meta))>1){
      stop("\"barcode\" and \"sampleID\" have to be included in arguement \"cell_meta\" as columns.")
    }


    # check if gene_of_interest is correct
    if(sum(!gene_of_interest%in%names(gene_list))==0){
      all_genes <- gene_of_interest
      gene_list <- gene_list[all_genes]
    }else{
      stop(" \"gene_of_interest\" is not valide. Please provide real gene names.")
    }



    ################### pre-step: reduce the genes according to the TSS clusters
    tss_clusters$start <- tss_clusters$org_start
    tss_clusters$end <- tss_clusters$org_end
    clusters_gr <- makeGRangesFromDataFrame(df = tss_clusters,
                                            keep.extra.columns = TRUE)
    resized_gene_list <- resize(gene_list,width(gene_list)+gene_extension,fix = "end")
    overlaps <- findOverlaps(resized_gene_list, clusters_gr)
    overlaps <- data.frame(overlaps)
    genes_idx <- unique(overlaps$queryHits)
    cluster_idx <- unique(overlaps$subjectHits)

    gene_list <- gene_list[genes_idx]
    tss_clusters <- tss_clusters[cluster_idx,]

    ## out reduced results: reduced genes
    all_genes <- names(gene_list)
    ## out reduced results：tss_clsuters
    clusters_gr <- makeGRangesFromDataFrame(df = tss_clusters,
                                            keep.extra.columns = TRUE)


    ######################## quantify: gene by gene
    print(paste0("Start quantifing ",length(clusters_gr)," clusters on ", length(gene_list)," genes"))
    tss_outs <- parallel:: mclapply(1:length(gene_list), function(iii){
      tryCatch({
        region <- gene_list[iii]
        temp_result <- query_based_on_single_regionInput(single_region = region,
                                                         query_regions = clusters_gr)
        outs <- TSS_quantify_perGene(region = region,
                                     temp_result = temp_result,
                                     reads_percentage = reads_percentage,
                                     read_number = read_number,
                                     gene_extension = gene_extension,
                                     bam_path = bam_path,
                                     cell_meta = cell_meta,
                                     if_paired = if_paired)
        if (iii%%500 == 0){
          gc()
        }

        return(outs)
      }, error=function(e){})

    },mc.cores = ncore)

    names(tss_outs) <- all_genes
    print("Counting reads finished!")
    tss_outs <- tss_outs[lengths(tss_outs)>0]
    tss_outs <- rbindlist(tss_outs,use.names=T)
    setDT(tss_outs,"TSS_clusters")

    setnames(tss_outs,
             cell_meta$barcode,
             paste0(cell_meta$sampleID,"_",cell_meta$barcode) )
    cell_meta$barcode <- paste0(cell_meta$sampleID,"_",cell_meta$barcode)

    return(list(TSS_clsuter_counts = tss_outs,
                Column_meta = cell_meta))
  }


}
