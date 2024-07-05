## a function to check the location on which exons

check_loc_exons <- function(exons, loc) {
  in_exon <- FALSE
  exon_idx <- NA
  after_which_exon <- NA


  # Check if loc is in any exon
  for (i in 1:nrow(exons)) {
    if (loc >= exons$start[i] && loc <= exons$end[i]) {
      in_exon <- TRUE
      exon_idx <- i
      break
    }
  }

  if(in_exon == FALSE){
    temp_dist <- exons$end - loc
    before_idx <- which(temp_dist < 0)
    after_which_exon <- which.max(before_idx)
  }



  outs <- data.frame(in_exon = in_exon,
                     exon_idx = exon_idx,
                     after_which_exon = after_which_exon)
  return(outs)
}
