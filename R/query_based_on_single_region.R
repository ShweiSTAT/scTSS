# a function to query from hit based on input single region

query_based_on_single_regionInput <- function(single_region,
                                              query_regions){
  overlaps <- suppressWarnings(findOverlaps(single_region, query_regions))
  overlaps <- data.frame(overlaps)
  mappedHit_idx <- unique(overlaps$subjectHits)
  query_regions_out <- query_regions[mappedHit_idx,]

  return(query_regions_out)
}
