
# to get the gene names for the genes with only one annotated TSS
getOneTSSgenes <- function(tx_list_by_gene, ncore){
  candiate_gene_names <- names(tx_list_by_gene)

  oneTSSgenes <- mclapply(candiate_gene_names,function(gene_names){
    tx_temp <- tx_list_by_gene[gene_names]
    gene_strand <- as.character(unique(strand(tx_temp)))
    if(gene_strand == "+"){
      tx_tss <- unique(unlist(start(tx_temp)))
    }else{
      tx_tss <- unique(unlist(end(tx_temp)))
    }
    ## record gene names for one-tss genes
    if(length(tx_tss) == 1){
      return(names(tx_temp))
    }
  },mc.cores=ncore)
  outs <- unlist(oneTSSgenes)

  return(outs)
  }


