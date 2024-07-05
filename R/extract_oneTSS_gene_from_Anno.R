# a function to extract oneTSS annotated TSS from .gtf file

extract_oneTSS_gene_from_Anno <- function(txdb_tobe_extracted){
  tx_list_by_gene <- transcriptsBy(txdb_tobe_extracted,by=c("gene"))

  # extract annotated promoters
  tss <- promoters(tx_list_by_gene, upstream = 0, downstream = 1)

  # remove duplicate promoter regions
  # tss_df <- data.frame(tss)
  # if("group_name" %in% colnames(tss_df)){
  #   unique_tss_df <- tss_df %>%
  #     distinct(group_name, start, end, .keep_all = TRUE)
  #   tss_count <- unique_tss_df %>%
  #     group_by(group_name) %>%
  #     summarize(tss_count = n())
  #   # extract genes with only one annotated TSSs
  #   single_tss_genes <- tss_count %>%
  #     filter(tss_count == 1) %>%
  #     select(group_name)
  #   single_tss_gene_names <- single_tss_genes$group_name
  # }else{
    # tss_loc_list <- start(tss)
    # tss_count <- lapply(tss_loc_list, function(x){
    #   length(unique(x))
    # })
    # tss_count <- unlist(tss_count)
    # single_tss_gene_names <- names(which(tss_count==1))
    # }

  tss_loc_list <- start(tss)
  tss_count <- lapply(tss_loc_list, function(x){
    length(unique(x))
  })
  tss_count <- unlist(tss_count)
  single_tss_gene_names <- names(which(tss_count==1))

  return(single_tss_gene_names)
}
