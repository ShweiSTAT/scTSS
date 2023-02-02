
## to read reads from .bam
reading_reads_loc <- function(region,bam_path,if_paired,gene_extension){
  if(if_paired){
    region<- resize(region,width(region)+gene_extension,fix = "end")
    param <- ScanBamParam(which = region,
                          tag = c("AN","CB"),
                          what = c("flag", "mrnm", "mpos"))
    reads <- readGAlignments(bam_path, param = param,use.names=TRUE)

    # drop AN reads and make pairs
    AN.info <- reads@elementMetadata$AN
    drop.idx.AN <- which( (!is.na(AN.info)) > 0 )
    if (length(drop.idx.AN)>0){
      reads <- reads[-drop.idx.AN]
    }
    readpairs <- makeGAlignmentPairs(reads,  strandMode=0, use.mcols=TRUE)
    readpairs <- data.frame(readpairs)

    # check if switch first and last needed
    switch.idx <- readpairs[,"start.last"] <= readpairs[,"start.first"]
    tp <- readpairs[switch.idx, c("start.last","end.last")]
    readpairs[switch.idx, c("start.last","end.last")] <- readpairs[switch.idx, c("start.first","end.first")]
    readpairs[switch.idx, c("start.first","end.first")] <- tp

    # drop reads that are out of the gene range
    readpairs <- readpairs[readpairs$start.first >= start(region)-gene_extension,]
    readpairs <- readpairs[readpairs$end.last <= end(region)+gene_extension,]

    return(readpairs)

  }else{
    # strandmode = 0
    # i = which(names(gene_list)==gene_id)
    # # read in reads
    # region <- gene_list[i]
    region<- resize(region,width(region)+gene_extension,fix = "end")
    param <- ScanBamParam(which = region,
                          tag = c("AN","CB"),
                          what = c("flag", "mrnm", "mpos"))
    reads <- readGAlignments(bam_path, param = param,use.names=TRUE)

    # drop AN reads
    AN.info <- reads@elementMetadata$AN
    drop.idx.AN <- which( (!is.na(AN.info)) > 0 )
    if (length(drop.idx.AN)>0){
      reads <- reads[-drop.idx.AN]
    }

    # drop reads that are out of the gene range
    reads <- data.frame(reads)
    reads <- reads[reads$start >= start(region)-gene_extension,]
    reads <- reads[reads$end <= end(region)+gene_extension,]
    return(reads)
  }
}
