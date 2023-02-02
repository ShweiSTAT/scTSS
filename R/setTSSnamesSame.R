## a function to set the rownames (TSS) the same across samples
## the input should be a list of usage matrices and separated by condotions
setTSSnamesSame <- function(x, join_type="inner"){

  if(sum(!(join_type=="inner"|join_type=="outter"))>0){
    stop("ERROR: Enter a correct joint_type" )
  }

  if(join_type == "inner"){
  allnames <- data.table(merged_tss=Reduce(intersect,lapply(x, function(y){y$merged_tss})))
  }else{
    allnames <- data.table(merged_tss=unique(unlist(lapply(x, function(y){y$merged_tss}))))
  }

  ## reset names
  x <- lapply(x, function(y){
    y <- merge(x=allnames,y=y,by="merged_tss",all.x=TRUE)
    removeNAs(y)
    setDT(y,key="merged_tss")
  })
  return(x)
}
