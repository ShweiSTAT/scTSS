## a function to perform DU test for each gene
DU_test <- function(test_gene,
                    a,
                    z,
                    z_per,
                    sample_label){
  # get tss candidates for this gene
  this_gene_tss <- unlist(test_gene[[1]][,1])
  # transferring data to a format that can be calculated as matrix
  test_gene <- lapply(test_gene, function(x){
    x <- x[,-1]
    tt <- colSums(x)
    x <- x[,-which(tt==0),with=FALSE]
    x <- data.table(x)
    return(x)
  })

  # a matrix to store DU information
  test_rst <- matrix(NA,nrow=length(this_gene_tss),ncol=7)
  rownames(test_rst) <- this_gene_tss
  colnames(test_rst) <- c("p_value","fold_change",
                          "max_label","max_meanUsage",
                          "min_label","min_meanUsage","num_of_cells")

  ## calculate the fold change between the lowest and highest in mean usage
  meanUsage <- lapply(test_gene, function(x){
    rowMeans(x)
  })
  meanUsage <- data.frame(do.call(rbind,meanUsage))
  colnames(meanUsage) <- this_gene_tss
  meanUsage$label <- sample_label$sample_label
  meanUsage <- aggregate(.~label, data=meanUsage, mean)
  for (TSS in this_gene_tss){
    tempLabel <- meanUsage$label
    tempMeanUsage <- meanUsage[,TSS]
    test_rst[TSS,2:7] <- c(max(tempMeanUsage)/min(tempMeanUsage),
                           tempLabel[which.max(tempMeanUsage)],
                           max(tempMeanUsage),
                           tempLabel[which.min(tempMeanUsage)],
                           min(tempMeanUsage),sum(lengths(test_gene)))

  }


  ## The major DU function
  H <- z%*%solve(t(z)%*%z)%*%t(z)
  # test on each tss
  for (tssIDX in 1:length(this_gene_tss)){
    # get a D matrix for every TSS
    D_list <- lapply(test_gene, function(x){
      return(x[tssIDX,])
    })
    # make a distance matrix
    nlist <- length(D_list)
    D <- array(0,c(nlist,nlist))
    for (i in 1:(nlist-1)){
      vali = as.numeric(as.character(D_list[[i]]))
      for (j in (i+1):nlist){
        valj = as.numeric(as.character(D_list[[j]]))
        theval <- round(wasserstein1d(vali, valj, p=1),digits = 10)
        D[i,j] <- D[j,i] <- theval
      }
    }

    # calculate other F_sudo
    A <- (-1/2)*D^2
    n <- nrow(z)
    centering <- (diag(n)-(1/n)*matrix(1,n,n))
    G <- centering%*%A%*%centering
    F_sudo <- sum(diag(H%*%G%*%H)) /(sum(diag( (diag(n)-H)%*%G%*%(diag(n)-H))))
    # permute Z and calculate the H0 distribution
    F_per <- numeric()
    for (ii in 1: a){
      Z_per <- z_per[[ii]]
      H_per <- Z_per%*%solve(t(Z_per)%*%Z_per)%*%t(Z_per)
      F_per <- c(F_per,
                 sum(diag(H_per%*%G%*%H_per)) /(sum(diag( (diag(n)-H_per)%*%G%*%(diag(n)-H_per)))) )
    }
    # calculate p_value
    p_val <- mean(F_sudo <= F_per)
    test_rst[tssIDX,1] <- p_val
  }
  return(test_rst)
}
