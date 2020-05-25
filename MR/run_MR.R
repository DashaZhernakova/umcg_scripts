library(TwoSampleMR)
library(MRInstruments)
library(MRPRESSO)


run_mr = function(exp_dat, out_dat, exp_table = NULL){
  if (is.null(out_dat)) return (NULL)
  
  dat <- harmonise_data(
    exposure_dat = exp_dat, 
    outcome_dat = out_dat
  )
  res <- mr(dat, method_list=c("mr_ivw", "mr_wald_ratio"))
  
  if (nrow(res) > 0){
    dat = dat[dat$mr_keep == T,]
    res$SNP <- paste(unique(dat$SNP), collapse = ',')
    if (!is.null(exp_table)){
      res$type <- paste(unique(exp_table[exp_table$SNP %in% dat$SNP,"type"]), collapse = ',')
    }
    res$egger_intercept_pval <- NA
    res$heterogeneity_Q_pval <- NA
    res$heterogeneity_Q <- NA
    #res$steiger_correct <- NA
    res$weighted_median_pval <- NA
    res$weighted_median_beta <- NA
    res$weighted_median_se <- NA
    res$egger_pval <- NA
    res$egger_beta <- NA
    res$egger_se <- NA
    res$mr_presso_pval <- NA
    res$mr_presso_outlier_cor_pval <- NA
    res$mr_presso_global <- NA
    res$leave_one_out_pval <- NA
    res$samplesize.outcome <-dat$samplesize.outcome[1]


    egger_intercept_pval <- mr_pleiotropy_test(dat)$pval
    if (!is.null(egger_intercept_pval)){
      res[, "egger_intercept_pval"] <- egger_intercept_pval
    }
    
    hetero<-  mr_heterogeneity(dat, method_list = "mr_ivw")
    if (nrow(hetero) > 0){
      res$heterogeneity_Q_pval <- hetero$Q_pval
      res$heterogeneity_Q <- hetero$Q
    }
    
    #steiger <- directionality_test(dat)$correct_causal_direction
    #if (!is.null(steiger)){  
    #  res[, "steiger_correct"] <- steiger
    #}
     
    
    pres = NULL
    if (res$nsnp > 3){
    pres <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)
    }
    if (!is.null(pres)){
        res$mr_presso_pval <- pres$`Main MR results`[1,6]
        res$mr_presso_outlier_cor_pval <- pres$`Main MR results`[2,6]
        res$mr_presso_global <- pres$`MR-PRESSO results`$`Global Test`$Pvalue
      }
    
    wm <- mr(dat, method_list=c("mr_weighted_median"))
    if (nrow(wm) > 0){  
      res$weighted_median_pval <- wm$pval
      res$weighted_median_beta <- wm$b
      res$weighted_median_se <- wm$se
    }
    
    egg <- mr(dat, method_list=c("mr_egger_regression"))
    if (nrow(egg) > 0){
      res$egger_pval <- egg$pval
      res$egger_beta <- egg$b
      res$egger_se <- egg$se
    }
    res$filter_loo = TRUE 
    res_loo <- mr_leaveoneout(dat)
    if (nrow(res_loo) > 1){
      res$leave_one_out_pval <- paste(res_loo[res_loo$SNP != "All","p"], collapse = ',')
    
    if (max(res_loo$p, na.rm = T) > 0.05){
      res$filter_loo=FALSE
    }
    }
    res$filter_all = TRUE
    if ((!is.na(res$heterogeneity_Q_pval) & res$heterogeneity_Q_pval < 0.05) | (!is.na(res$egger_intercept_pval) & res$egger_intercept_pval < 0.05) | (!is.na(res$weighted_median_pval) & res$weighted_median_pval > 0.05) | (!is.na(res$mr_presso_pval) & res$mr_presso_pval > 0.05) | (!is.na(res$mr_presso_global) & res$mr_presso_global < 0.05 & res$mr_presso_outlier_cor_pval > 0.05)){
      res$filter_all = FALSE
    }
  }
  return(res)
}
