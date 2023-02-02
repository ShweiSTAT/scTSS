## a function to remove NAs for data.table
removeNAs <- function(DT) {
  for (j in names(DT))
    set(DT,which(is.na(DT[[j]])),j,0)
  for (j in seq_len(ncol(DT)))
    set(DT,which(is.na(DT[[j]])),j,0)
}
