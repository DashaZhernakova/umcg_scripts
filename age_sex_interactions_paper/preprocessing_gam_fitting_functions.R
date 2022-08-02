library(RColorBrewer)
library('dplyr')
library('mgcv')
#library('gamlss')

# Combine phenotype and trait of interest and remove outliers in the trai of interest
rm_na_outliers <- function(traits_m, pheno_m, idx, method = "zscore", scale_tr = F, log_tr = F, int_tr = F){
  traits_na.rm <- traits_m[!is.na(traits_m[,idx]),idx]
  pheno_na.rm <- pheno_m[!is.na(traits_m[,idx]),]
  all(row.names(traits_na.rm) == row.names(pheno_na.rm))
  merged_tab <- cbind(traits_na.rm, pheno_na.rm)
  colnames(merged_tab) <- c(colnames(traits_m)[idx], colnames(pheno_na.rm))
  row.names(merged_tab) <- row.names(pheno_na.rm)
  
  pheno_is_factor <- F
  if (length(unique(merged_tab[,1])) < 3) pheno_is_factor <- T
  
  # apply transfromations
  if (pheno_is_factor & any(scale_tr, log_tr, int_tr)){
    scale_tr = F
    log_tr = F
    int_tr = F
    pheno_is_factor <- T
    message("Outcome is a factor, skipping all transformations!")
  }
  
  #min_val <- min(merged_tab[merged_tab[,1] != 0,1], na.rm = T)
  min_val = 1
  if (log_tr & scale_tr){
    merged_tab[,1] <- scale(log(merged_tab[,1] + min_val))
  } else if (log_tr & ! scale_tr) {
    merged_tab[,1] <- log(merged_tab[,1] + min_val)
  } else if (!log_tr & scale_tr) {
    merged_tab[,1] <- scale(merged_tab[,1])
  } else if (int_tr) {
    merged_tab[,1] <- qnorm((rank(merged_tab[,1],na.last="keep")-0.5)/sum(!is.na(merged_tab[,1])))
  }
  
  
  # remove outliers ( skip this if outcome is a factor )
  if (! pheno_is_factor){
    
    if (method == "zscore"){ 
      # Zscore < 3
      #tab_nooutliers <- rbind(w[abs(w[,1] - mean(w[,1]))/sd(w[,1]) < 3,], m[abs(m[,1] - mean(m[,1]))/sd(m[,1]) < 3,])
      tab_nooutliers <- merged_tab[abs(merged_tab[,1] - mean(merged_tab[,1])) / sd(merged_tab[,1]) < 3,]
    } else if (method == "IQR"){
      # less than 1.5 IQR from the 1st and 3rd quantiles
      q1 <- quantile(merged_tab[,1], probs = 0.25)
      q3 <- quantile(merged_tab[,1], probs = 0.75)
      iqr <- q3 - q1
      tab_nooutliers <- merged_tab[(merged_tab[,1] < (q3 + 3*iqr)) & (merged_tab[,1] > (q1 - 3*iqr)),]
      
    } else {
      message ("Wrong method! No outlier removal!")
      tab_nooutliers <- merged_tab
    }
    
  } else {
    tab_nooutliers <- merged_tab
  }
  
  
  num_outliers <- (nrow(merged_tab) - nrow(tab_nooutliers))
  perc_outliers <- format(100*num_outliers/nrow(merged_tab), digits = 3)
  cat("\tRemoved", num_outliers, "(", perc_outliers, "%)","outliers\n")
  return(tab_nooutliers)
}

# Fit a GAM with age : gender interaction and (optional) correction for covariates (linear or spline)
plot_scatter_and_gam2 <- function(merged_tab, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots, ymax_hist = 1, label = '', gam_family = gaussian(), min_age = 20, max_age = 80, add_inter_p_to_plot = T, plot_title = NULL, interp_cutoff = 1, plot_points = T, add_breakpoints = F, t_threshold = 3, derivatives_cutoff = 0.0002, ylims = NULL, log_tr = F, plot_density = F){  
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
  
  if (length(unique(merged_tab[,1])) < 3 ){
     cat("Outcome is binary, changing family to binomial with logit\n")
     family=binomial(link = 'logit')
  }
  res_dif = NULL
  if (length(covariates_linear) == 0 & length(covariates_nonlinear) == 0){
    
    gam.fit <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML", family = gam_family)
    gam.fit.ordered_sex <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
                 data = merged_tab, method = 'REML', family = gam_family)
    gam.fit.ordered_sex.no_interaction <- gam(phenotype ~ ord_gender_F1M2 + s(age), 
                               data = merged_tab, method = 'REML', family = gam_family)
    
    gam.p <- summary(gam.fit.ordered_sex)$s.pv[length(summary(gam.fit.ordered_sex)$s.pv)]
    gam.g_beta <- gam.fit.ordered_sex$coefficients["ord_gender_F1M2.L"]
    gam.g_pv <- summary(gam.fit.ordered_sex)$p.pv["ord_gender_F1M2.L"]
    gam.cohen_f2 <- calculate_cohens_f2(gam.fit.ordered_sex, gam.fit.ordered_sex.no_interaction)
    
    #if (gam.p > interp_cutoff){
    #  return (list("dif" = NULL, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv, "cohen_f2" = gam.cohen_f2))
    #}
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points), gender_F1M2 = c('1', '2'))) 
    new.y <- data.frame(predict(gam.fit, newdata = new.x, se.fit = TRUE, type = "response"))
    pdat <- data.frame(new.x, new.y)
    pdat <- rename(pdat, pred = fit, SE = se.fit)
    pdat <- mutate(pdat, lwr = pred - 1.96 * SE, upr = pred + 1.96 * SE) # calculating the 95% confidence interval
    
    
  } else { # Correct for covariates
    terms_linear_covar <- ""
    terms_nonlinear_covar <- ""
    if (length(covariates_linear) > 0) terms_linear_covar <- paste0("+", paste(covariates_linear, collapse = "+"))
    if (length(covariates_nonlinear) > 0) terms_nonlinear_covar <- paste0("+ s(", paste(covariates_nonlinear, collapse = ")+ s("), ")")
    
    gam.fit <- gam(as.formula(paste("phenotype ~ gender_F1M2 ", terms_linear_covar, terms_nonlinear_covar,
                               "+ s(age) + s(age, by = gender_F1M2)", sep = " ")), 
              data = merged_tab, method = "REML")
    gam.fit.ordered_sex <- gam(as.formula(paste("phenotype ~ ord_gender_F1M2 ", terms_linear_covar, terms_nonlinear_covar,
                               "+ s(age) + s(age, by = ord_gender_F1M2)", sep = " ")), 
              data = merged_tab, method = "REML")
    gam.fit.ordered_sex.no_interaction <- gam(as.formula(paste("phenotype ~ ord_gender_F1M2 ", terms_linear_covar, terms_nonlinear_covar,
                                                               "+ s(age)", sep = " ")), 
                                              data = merged_tab, method = "REML")
    
    gam.p <- summary(gam.fit.ordered_sex)$s.pv[length(summary(gam.fit.ordered_sex)$s.pv)]
    gam.g_beta <- gam.fit.ordered_sex$coefficients["ord_gender_F1M2.L"]
    gam.g_pv <- summary(gam.fit.ordered_sex)$p.pv["ord_gender_F1M2.L"]
    gam.cohen_f2 <- calculate_cohens_f2(gam.fit.ordered_sex, gam.fit.ordered_sex.no_interaction)
    
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points), gender_F1M2 = c('1', '2'))) 
    for (c in c(covariates_linear, covariates_nonlinear)){
      if (is.factor(merged_tab[,c])){
        new.x[,c] <- 0
        new.x[,c] <- as.factor(new.x[,c])
      } else {
        new.x[,c] <- mean(merged_tab[,c])
      }
    }
    
    new.y <- data.frame(predict(gam.fit, newdata = new.x, se.fit = TRUE, type = "response"))
    pdat <- data.frame(new.x, new.y)
    pdat <- rename(pdat, pred = fit, SE = se.fit)
    pdat <- mutate(pdat, lwr = pred - 1.96 * SE, upr = pred + 1.96 * SE) # calculating the 95% confidence interval
  }
  
  if (log_tr){
    merged_tab[,1] <-  exp(merged_tab[,1]) - 1
    pdat$pred <- exp(pdat$pred) - 1
    pdat$lwr <- exp(pdat$lwr) - 1
    pdat$upr <- exp(pdat$upr) - 1
  }

  #if (gam.p > interp_cutoff){
  #  return (list("dif" = NULL, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv, "cohen_f2" = gam.cohen_f2))
  #}
  p <- NULL 
  if (make_plots & gam.p < interp_cutoff & plot_density == F){
     draw_plot(merged_tab, pheno_name, pdat, gam.p, min_age = min_age, max_age = max_age, add_inter_p_to_plot = add_inter_p_to_plot, plot_title = plot_title, plot_points = plot_points,  ymax_hist = ymax_hist, label = paste0("Cohen's f^2 = ", formatC(gam.cohen_f2, digits = 4, format = "f")), ylims = ylims)
     #draw_smooth_scatter(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot = add_inter_p_to_plot, plot_title = plot_title, label = paste0("Cohen's f^2 = ", formatC(gam.cohen_f2, digits = 4, format = "f")))
  } else if (make_plots & gam.p < interp_cutoff & plot_density){
    p <- draw_contour_plot(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot, plot_title,  ymax_hist, label = paste0("Cohen's f^2 = ", formatC(gam.cohen_f2, digits = 4, format = "f")), ylims = ylims)
  }
  
  return (list("pdat" = pdat, "dif" = res_dif$diff, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv, "cohen_f2" = gam.cohen_f2, "p" = p))
} 




# Fit a GAM with age : gender interaction and (optional) correction for covariates (linear or spline)
fit_gam_age_only <- function(merged_tab, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots = T,  gam_family = gaussian(), min_age = 20, max_age = 80, log_tr = F){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  if (length(unique(merged_tab[,1])) < 3 ){
     cat("Outcome is binary, changing family to binomial with logit\n")
     family=binomial(link = 'logit')
  }
 
  if (length(covariates_linear) == 0 & length(covariates_nonlinear) == 0){
    gam.fit <- gam(phenotype ~ s(age), data = merged_tab, method = "REML", family = gam_family)
    
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points))) 
    new.y <- data.frame(predict(gam.fit, newdata = new.x, se.fit = TRUE, type = "response"))
    pdat <- data.frame(new.x, new.y)
    pdat <- rename(pdat, pred = fit, SE = se.fit)
    pdat <- mutate(pdat, lwr = pred - 1.96 * SE, upr = pred + 1.96 * SE) # calculating the 95% confidence interval
    
    
  } else { # Correct for covariates
    terms_linear_covar <- ""
    terms_nonlinear_covar <- ""
    if (length(covariates_linear) > 0) terms_linear_covar <- paste0("+", paste(covariates_linear, collapse = "+"))
    if (length(covariates_nonlinear) > 0) terms_nonlinear_covar <- paste0("+ s(", paste(covariates_nonlinear, collapse = ")+ s("), ")")
    
    gam.fit <- gam(as.formula(paste("phenotype ~ ", terms_linear_covar, terms_nonlinear_covar,
                                    "+ s(age)", sep = " ")), 
                   data = merged_tab, method = "REML")
    
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points))) 
    for (c in c(covariates_linear, covariates_nonlinear)){
      if (is.factor(merged_tab[,c])){
        new.x[,c] <- 0
        new.x[,c] <- as.factor(new.x[,c])
      } else {
        new.x[,c] <- mean(merged_tab[,c])
      }
    }
    
    new.y <- data.frame(predict(gam.fit, newdata = new.x, se.fit = TRUE, type = "response"))
    pdat <- data.frame(new.x, new.y)
    pdat <- rename(pdat, pred = fit, SE = se.fit)
    pdat <- mutate(pdat, lwr = pred - 1.96 * SE, upr = pred + 1.96 * SE) # calculating the 95% confidence interval
  }
    if (log_tr){
    pdat$pred <- exp(pdat$pred) - 1
    pdat$lwr <- exp(pdat$lwr) - 1
    pdat$upr <- exp(pdat$upr) - 1
  }

  return (pdat)
}


# Fit a GAM with age : gender interaction and (optional) correction for covariates (linear or spline)
fit_gam_age_gender_only <- function(merged_tab, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots = T,  gam_family = gaussian(), min_age = 20, max_age = 80, log_tr = F){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  if (length(unique(merged_tab[,1])) < 3 ){
    cat("Outcome is binary, changing family to binomial with logit\n")
    family=binomial(link = 'logit')
  }
  merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
  
  if (length(covariates_linear) == 0 & length(covariates_nonlinear) == 0){
    gam.fit.both <- gam(phenotype ~ gender_F1M2 + s(age), data = merged_tab, method = "REML", family = gam_family)
    gam.fit.age <- gam(phenotype ~ s(age), data = merged_tab, method = "REML", family = gam_family)
    gam.fit.sex<- gam(phenotype ~ gender_F1M2 , data = merged_tab, method = "REML", family = gam_family)  
  } else { # Correct for covariates
    terms_linear_covar <- ""
    terms_nonlinear_covar <- ""
    if (length(covariates_linear) > 0) terms_linear_covar <- paste0("+", paste(covariates_linear, collapse = "+"))
    if (length(covariates_nonlinear) > 0) terms_nonlinear_covar <- paste0("+ s(", paste(covariates_nonlinear, collapse = ")+ s("), ")")
    
    gam.fit.both <- gam(as.formula(paste("phenotype ~ ", terms_linear_covar, terms_nonlinear_covar,
                                    "+ gender_F1M2 + s(age)", sep = " ")), 
                   data = merged_tab, method = "REML")
    gam.fit.age <- gam(as.formula(paste("phenotype ~ ", terms_linear_covar, terms_nonlinear_covar,
                                         " + s(age)", sep = " ")), 
                        data = merged_tab, method = "REML")
    gam.fit.sex <- gam(as.formula(paste("phenotype ~ ", terms_linear_covar, terms_nonlinear_covar,
                                         "+ gender_F1M2", sep = " ")), 
                        data = merged_tab, method = "REML")
  }
  p_age <- summary(gam.fit.both)$s.table["s(age)",4]
  f2_age <- calculate_cohens_f2(gam.fit.both, gam.fit.sex)
  
  p_sex <- summary(gam.fit.both)$p.table["gender_F1M22",4]
  f2_sex <- calculate_cohens_f2(gam.fit.both, gam.fit.age)
  return (c(p_age, f2_age, p_sex, f2_sex))
}

calc_explained_variance_CV <- function(merged_tab, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), gam_family = gaussian(), min_age = 20, max_age = 80){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  if (length(unique(merged_tab[,1])) < 3 ){
    cat("Outcome is binary, changing family to binomial with logit\n")
    family=binomial(link = 'logit')
  }
  merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
  
  n = nrow(merged_tab)
  k = 5
  rep_n <- 10
  res_table <- data.frame(matrix(nrow = k*rep_n, ncol = 4))
  
  cnt <- 1
  for (j in 1:rep_n){
    shuffled <- merged_tab[sample(n),]
    folds <- cut(seq(1,n),breaks=k,labels=FALSE)
    
    for(i in 1:k){
      testInd <- which(folds==i,arr.ind=TRUE)
      test_data <- shuffled[testInd, ]
      train_data <- shuffled[-testInd, ]
      
      res_table[cnt,] <- calc_explained_variance_from_gam(train_data, test_data, pheno_name, covariates_linear, covariates_nonlinear, gam_family)
      cnt <- cnt + 1
    }
  }
  mean_res <- colMeans(res_table)
  exp_var_dif <- c(mean_res[1], mean_res[2:4] - mean_res[1:3])
  
  return(exp_var_dif)
}

