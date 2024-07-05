# a function to conduct glmm DU test for each gene

glmm_test_for_a_gene <- function(temp_full_data,
                                 col_meta,
                                 agg_level = "bulk",
                                 full_model,
                                 base_model){

  gene_counts <- colSums(temp_full_data)

  ## prepare for design matrix
  datx <- data.frame(col_meta,
                     gene_counts = gene_counts)

  pvals <- rep(NA,nrow(temp_full_data))
  betas <- list()
  max_usage_con <- rep(NA,nrow(temp_full_data))
  max_usage <- rep(NA,nrow(temp_full_data))
  min_usage_con <- rep(NA,nrow(temp_full_data))
  min_usage <- rep(NA,nrow(temp_full_data))
  gene_exp_level <- rep(NA,nrow(temp_full_data))
  TSS_exp_level <- rep(NA,nrow(temp_full_data))
  for(TSSidx in 1:nrow(temp_full_data)){

    tryCatch({
      temp_TSS <- temp_full_data[TSSidx,]
      dat <- data.frame(TSS_counts = temp_TSS,
                        datx)
      dat <- dat[which(dat$gene_counts > 0),]
      ################
      ## out stats
      gene_exp_level[TSSidx] <- nrow(dat)
      TSS_exp_level[TSSidx] <- sum(dat$TSS_counts>0)

      condition_bentch <-  unique(datx$condition)
      con_usage <- rep(NA,length(condition_bentch))
      names(con_usage) <- condition_bentch
      for(con_idx in 1:length(condition_bentch)){
        temp_con <- condition_bentch[con_idx]
        temp_con_dat <- dat[which(dat$condition == temp_con),]
        con_usage[con_idx] <- sum(temp_con_dat$TSS_counts)/sum(temp_con_dat$gene_counts)
      }

      max_usage_con[TSSidx] <- condition_bentch[which.max(con_usage)]
      max_usage[TSSidx] <- max(con_usage)
      min_usage_con[TSSidx] <- condition_bentch[which.min(con_usage)]
      min_usage[TSSidx] <- min(con_usage)

      #################

      if(agg_level =="bulk"){
        dat.df <- stats::aggregate(cbind(TSS_counts,gene_counts) ~ sampleID+condition, dat, sum)
      }else{
        dat.df <- dat
      }
      dat <- NULL

      suppressWarnings(mod.full <- glmer(formula = full_model,
                        data=dat.df,
                        verbose = FALSE,
                        family = "binomial"))
      suppressWarnings(mod.base <- glmer(formula = base_model,
                        data=dat.df,
                        verbose = FALSE,
                        family = "binomial"))

      rst <- stats::anova(mod.full, mod.base, test ="Chisq")
      pvals[TSSidx] <- rst$`Pr(>Chisq)`[2]
      betas[[TSSidx]] <- coef(summary(mod.full))[,1]
    }, error=function(e){})

  }
  outs <- data.frame(TSS_names = rownames(temp_full_data),
                     gene_exp_level = gene_exp_level,
                     TSS_exp_level = TSS_exp_level,

                     pval = pvals,
                     #do.call(rbind,betas),
                     max_usage_con = max_usage_con,
                     max_usage = max_usage,
                     min_usage_con = min_usage_con,
                     min_usage = min_usage)


  if(length(grep("ntercept",colnames(outs))) ==1){
    colnames(outs)[grep("ntercept",colnames(outs))] <- "Intercept"
  }


  return(outs)
}
