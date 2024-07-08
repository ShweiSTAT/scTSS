# A function to split the original Near-site genomic regions
# for merged TSS clusters base on singe-end data

merge_match_for_SingleEnd_data <- function(single_unmerged_region,merged_regions){
  matched_merged_regions <- query_based_on_single_regionInput(single_region = single_unmerged_region,
                                                              query_regions = merged_regions)

  ## get gene info
  strand <- as.character(strand(single_unmerged_region))

  ## split regions
  matched_merged_regions <- as.data.frame(matched_merged_regions)
  single_unmerged_region <- as.data.frame(single_unmerged_region)

  if(strand == "-"){
    matched_merged_regions <- matched_merged_regions[order(matched_merged_regions$tss5prime,
                                                           decreasing = TRUE),]
  }
  if(nrow(matched_merged_regions)<=1){
    # if only one match, output the single unmerged regions
    outs <- data.frame(seqnames = single_unmerged_region$seqnames,
                       tss5prime = single_unmerged_region$tss5prime,
                       width = single_unmerged_region$width,
                       strand = single_unmerged_region$strand,
                       org_start = single_unmerged_region$org_start,
                       org_end = single_unmerged_region$org_end)
  }else{
    # if there are more than 1 match, split the orginal regions (before adj)
    split_rst_org_regions <- split_region(chr = single_unmerged_region$seqnames,
                                          start = single_unmerged_region$org_start,
                                          end = single_unmerged_region$org_end,
                                          strand = single_unmerged_region$strand,
                                          widths = matched_merged_regions$width)
    split_rst_org_regions <- data.frame(split_rst_org_regions)
    outs <- data.frame(seqnames = single_unmerged_region$seqnames,
                       tss5prime = matched_merged_regions$tss5prime,
                       width = matched_merged_regions$width,
                       strand = single_unmerged_region$strand,
                       org_start = split_rst_org_regions$start,
                       org_end = split_rst_org_regions$end)
    outs <- outs[order(outs$tss5prime,
                       decreasing = FALSE),]
  }

  return(outs)

}
