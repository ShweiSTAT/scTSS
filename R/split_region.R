## a function to split genomic regions from the 5' based on
## the intended width for each split regions

split_region <- function(chr, start, end, strand, widths) {
  if (strand == "+") {
    starts <- start + cumsum(c(0, widths[-length(widths)]))
    ends <- starts + widths - 1
  } else if (strand == "-") {
    ends <- end - cumsum(c(0, widths[-length(widths)]))
    starts <- ends - widths + 1
  } else {
    stop("Strand must be either '+' or '-'")
  }

  # Create GRanges object
  gr <- GRanges(seqnames = chr,
                ranges = IRanges(start = starts, end = ends),
                strand = strand)
  return(gr)
}
