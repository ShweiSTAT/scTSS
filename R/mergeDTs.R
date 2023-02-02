## a function to merge data.tables
mergeDTs <- function(dt_list, by = NULL) {
  Reduce(
    function(...) {
      merge(..., by = by, all = TRUE)
    }, dt_list)
}
