# a function to get mapped TSS table (gene by gene)

getMappedTSS_perGene <- function(region_of_this_gene,
                               bam_path,
                               if_paired_this_data,
                               gene_extension = 500){

  chr <- seqnames(region_of_this_gene)
  strand <- as.character(strand(region_of_this_gene))
  start <- start(region_of_this_gene)
  end <- end(region_of_this_gene)

  ## read in all the reads for a gene
  read_in <- reading_reads_loc(region = region_of_this_gene,
                               bam_path = bam_path ,
                               if_paired = if_paired_this_data,
                               gene_extension = gene_extension)

  ## tabulate the Mapped TSS
  if(if_paired_this_data){
    if (strand == "+"){
      temp <- table(read_in$start.first)
      MappedTSS_tab <- data.frame(seqnames = as.character(chr),
                             tss = as.numeric(names(temp)),
                             strand = strand,
                             tags = as.numeric(temp)
      )

    }else{
      temp <- table(read_in$end.last)
      MappedTSS_tab <- data.frame(seqnames = chr,
                             tss = as.integer(names(temp)),
                             strand = strand,
                             tags = as.numeric(temp)
      )

    }

  }else{
    if (strand == "+"){
      temp <- table(read_in$start)
      MappedTSS_tab <- data.frame(seqnames = as.character(chr),
                                  tss = as.numeric(names(temp)),
                                  strand = strand,
                                  counts = as.numeric(temp)
      )

    }else{
      temp <- table(read_in$end)
      MappedTSS_tab <- data.frame(seqnames = chr,
                                  tss = as.integer(names(temp)),
                                  strand = strand,
                                  counts = as.numeric(temp)
      )

    }
  }

  return(MappedTSS_tab)

}
