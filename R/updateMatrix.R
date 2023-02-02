## a function to update unmerged TSS count matrices to
## merged TSS matrices by dictionary
updateMatrix <- function(mat,merging_dictionary) {
  setDT(mat,key="unmerged_tss")
  mat <- merge(x=merging_dictionary,y=mat,key="unmerged_tss",all.y=T,sort=T)
  setDT(mat,key="merged_tss")
  removeNAs(mat)
  mat$unmerged_tss <- NULL
  return(mat)
}
