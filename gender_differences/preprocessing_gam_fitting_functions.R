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
  
  min_val <- min(merged_tab[merged_tab[,1] != 0,1], na.rm = T)
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
    w <- merged_tab[merged_tab$gender_F1M2 == 1,]
    m <- merged_tab[merged_tab$gender_F1M2 == 2,]
   
    if (method == "zscore"){ 
      # Zscore < 3
      tab_nooutliers <- rbind(w[abs(w[,1] - mean(w[,1]))/sd(w[,1]) < 3,], m[abs(m[,1] - mean(m[,1]))/sd(m[,1]) < 3,])
    } else if (method == "IQR"){
      # less than 1.5 IQR from the 1st and 3rd quantiles
      mq1 <- quantile(m[,1], probs = 0.25)
      mq3 <- quantile(m[,1], probs = 0.75)
      miqr <- mq3 - mq1
      m_clean <- m[m[,1] < mq3 + 1.5*miqr & m[,1] > mq1 - 1.5*miqr,]
      
      wq1 <- quantile(w[,1], probs = 0.25)
      wq3 <- quantile(w[,1], probs = 0.75)
      wiqr <- wq3 - wq1
      w_clean <- w[w[,1] < wq3 + 1.5*wiqr & w[,1] > wq1 - 1.5*wiqr,]
      
      tab_nooutliers <- rbind(w_clean, m_clean)
    } else {
      message ("Wrong method! No outlier removal!")
      tab_nooutliers <- merged_tab
    }
    
  } else {
    tab_nooutliers <- merged_tab
  }
  
  # convert binary phenotypes to factors
  #factor_cols <- which(sapply(tab_nooutliers, function(col) length(unique(col)) < 3))
  #for (c in factor_cols){
  #  if (c != 1){
  #  tab_nooutliers[,c] <- as.factor(tab_nooutliers[,c])
  #  }
  #}
  return(tab_nooutliers)
}

# Fit a GAM with age : gender interaction and (optional) correction for covariates (linear or spline)
plot_scatter_and_gam2 <- function(merged_tab, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots, label = '', gam_family = gaussian(), min_age = 20, max_age = 80, add_inter_p_to_plot = T, plot_title = NULL, interp_cutoff = 1, plot_points = T, add_breakpoints = F, t_threshold = 3, derivatives_cutoff = 0.0002){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
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
    
    if (gam.p > interp_cutoff){
      return (list("dif" = NULL, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv, "cohen_f2" = gam.cohen_f2))
    }
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
  
  if (gam.p > interp_cutoff){
    return (list("dif" = NULL, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv, "cohen_f2" = gam.cohen_f2))
  }
  breakpoints_intervals <- NULL
  breakpoints <- NULL
  if (add_breakpoints){
      breakpoints_intervals <- get_breakpoints(merged_tab, correctForCellCounts, t_threshold = t_threshold, derivatives_cutoff = derivatives_cutoff)
      #breakpoints2 <- get_breakpoints_derivatives(merged_tab, correctForCellCounts)
  }
 
  
  if (make_plots & gam.p < interp_cutoff){
    draw_plot(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot, plot_title, plot_points, breakpoints, breakpoints_intervals, label = paste0("Fsq = ", formatC(gam.cohen_f2, digits = 4, format = "f")))
  }
  
  return (list("pdat" = pdat, "dif" = res_dif$diff, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv, "breakpoints_intervals" = breakpoints_intervals, "cohen_f2" = gam.cohen_f2))
} 

# Calculate a simple difference between every fitted value pair in men and women
simple_diff <- function(pdat){
  pdat <- pdat[order(pdat$gender_F1M2, pdat$age),]
  dif <- pdat[pdat$gender_F1M2 == 2, "pred"] - pdat[pdat$gender_F1M2 == 1, "pred"]
  res_dif <- data.frame(age = pdat[pdat$gender_F1M2 == 1, "age"], diff = dif)
  return(res_dif)
}


# Calculate the difference between men and women with confidence intervals
# source: https://www.fromthebottomoftheheap.net/2017/10/10/difference-splines-i/
smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  res_dif <- data.frame(pair = paste(f1, f2, sep = '-'),
                        diff = dif,
                        se = se,
                        upper = upr,
                        lower = lwr)
}

plot_scatter_and_gamlss <- function(merged_tab, pheno_name, n_points = 300, make_plots = T, gam_family = NO(), label = '', min_age = 18, max_age = 75){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  
  res <- test_polynomial_interaction(merged_tab, gam_family)  
  inter_pval <- res[["lr_p.val"]]
  gam.fit <- res[["mod"]] 
  if (inter_pval > 0.05){
      return (NULL)
  }

  pdat <- with(merged_tab, expand.grid(age = seq(min_age, 75, length = n_points), 
                                       gender_F1M2 = c('1', '2')))
  pdat <- transform(pdat, pred = predict(gam.fit, newdata = pdat, type = "response"))
  pdat_dif <- expand.grid(age = seq(min_age, max_age, length = n_points),
                          gender_F1M2 = c('1', '2'))
  res_dif <- simple_diff(pdat)
  
  ylims = c(min(merged_tab[,1], pdat$pred), max(merged_tab[,1], pdat$pred))
  if (make_plots){

    ylabel <- pheno_name
    if (nchar(pheno_name) > 40){
      spl <- unlist(strsplit(pheno_name, "\\|"))
      ylabel <- spl[length(spl)]
      cex_main <- 0.8
      if (nchar(pheno_name) > 50) cex_main <- 0.7
      if (nchar(pheno_name) > 60) cex_main <- 0.6
      
    }
    
    
    ## draw base plot
    palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
    plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
         main = paste0(pheno_name, ' ', label, "\nLR.test interaction p = ", format(inter_pval, digits = 3)), 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, ylim = ylims)
    
    levs <- levels(merged_tab$gender_F1M2)
    cols = c("indianred1", "dodgerblue1")
    
    ## add the fitted lines
    for (l in seq_along(levs)) {
      dd <- pdat[pdat$gender_F1M2 == levs[l],]
      lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
      
    }
    
    #Plot the difference
    plot(res_dif$age, res_dif$diff, type = 'l')
  }
  
}

# Fit a GAM with age : gender interaction and (optional) correction for cell counts
ggplot_scatter_and_gam2 <- function(merged_tab, pheno_name, correctForCellCounts, n_points = 300, make_plots, label = '', gam_family = gaussian(), min_age = 18, max_age = 75){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  #pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  
  women <- merged_tab$gender_F1M2 == 1
  res_dif = NULL
  
  if (!correctForCellCounts){
    
    gam.fit <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML", family = gam_family)
    gam.fit.ordered_sex <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
               data = merged_tab, method = 'REML', family = gam_family)
    gam.nointer <- gam(phenotype ~ ord_gender_F1M2 + s(age), 
               data = merged_tab, method = 'REML', family = gam_family)
    
    gam.p <- summary(gam.fit.ordered_sex)$s.pv[length(summary(gam.fit.ordered_sex)$s.pv)]
    gam.g_beta <- gam.fit.ordered_sex$coefficients["ord_gender_F1M2.L"]
    gam.g_pv <- summary(gam.fit.ordered_sex)$p.pv["ord_gender_F1M2.L"]
    calculate_cohens_f2
    if (gam.p > 0.05){
      return (list("dif" = NULL, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv))
    }
    pdat <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points), 
                                         gender_F1M2 = c('1', '2')))
    pdat <- transform(pdat, pred = predict(gam.fit, newdata = pdat, type = "response"))
    
    
    
    if (! "Ordered Categorical" %in% gam_family$family){
      pdat_dif <- expand.grid(age = seq(min_age, max_age, length = n_points),
                              gender_F1M2 = c('1', '2'))
      #res_dif <- smooth_diff(gam.fit, pdat_dif, '2', '1', "gender_F1M2")
      #res_dif$age = seq(min_age, max_age, length = n_points)
      res_dif <- simple_diff(pdat)
    } else { # stupid way to go
      gam.fit <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
      pdat <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points), 
                                           gender_F1M2 = c('1', '2')))
      pdat <- transform(pdat, pred = predict(gam.fit, newdata = pdat, type = "response"))
      
    }
  } else { # Correct for cell counts
    
    gam.fit <- gam(phenotype ~ gender_F1M2 + ba + eo + er + gr + 
                ly + mo + tr + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
    
    gam.fit.ordered_sex <- gam(phenotype ~ ord_gender_F1M2 + s(age) + ba + eo + er + gr + 
                 ly + mo + tr + s(age, by = ord_gender_F1M2), 
               data = merged_tab, method = 'REML')
    gam.p <- summary(gam.fit.ordered_sex)$s.pv[length(summary(gam.fit.ordered_sex)$s.pv)]
    gam.g_beta <- gam.fit.ordered_sex$coefficients["ord_gender_F1M2.L"]
    gam.g_pv <- summary(gam.fit.ordered_sex)$p.pv["ord_gender_F1M2.L"]
    
    if (gam.p > 0.05){
      return (list("dif" = NULL, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv))
    }
    
    pdat <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points), gender_F1M2 = c('1', '2'), 
                                         ba = mean(ba), eo = mean(eo), er = mean(er), gr = mean(gr), 
                                         ly = mean(ly),  mo = mean(mo), tr = mean(tr)))
    pdat <- transform(pdat, pred = predict(gam.fit, newdata = pdat, type = "response"))
    
    if (! "Ordered Categorical" %in% gam_family$family){
      pdat_dif <- expand.grid(age = seq(min_age, max_age, length = n_points), gender_F1M2 = c('1', '2'), 
                              ba = mean(merged_tab$ba), eo = mean(merged_tab$eo), er = mean(merged_tab$er), gr = mean(merged_tab$gr), 
                              ly = mean(merged_tab$ly),  mo = mean(merged_tab$mo), tr = mean(merged_tab$tr))
      #res_dif <- smooth_diff(gam.fit, pdat_dif, '2', '1', "gender_F1M2")
      #res_dif$age = seq(min_age, max_age, length = n_points)
      res_dif <- simple_diff(pdat)
    } else { # stupid way to go
      gam.fit <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
      pdat <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points), 
                                           gender_F1M2 = c('1', '2')))
      pdat <- transform(pdat, pred = predict(gam.fit, newdata = pdat, type = "response"))
    }
  }
  
  
  ylims <- with(merged_tab, range(phenotype))
  ## draw base plot
  palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
  
  p1 <- ggplot(merged_tab, aes(y = phenotype, x = age,  col = gender_F1M2)) + geom_point(aes(col = gender_F1M2)) + 
    scale_color_manual(values = c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125),"indianred1", "dodgerblue1")) + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
          axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid')) + 
    ylab(pheno_name) + 
    ggtitle(paste0(pheno_name, ' ', label, "\nGAM interaction p = ", format(gam.p, digits = 3))) + 
    scale_x_discrete(name ="age", limits=seq(10,80,10))
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    p1 <- p1 + geom_line(data = dd, aes(x= age, y = pred), colour = cols[[l]], size = 1)
  }
  
  #Plot the difference
  if (! is.null(res_dif = NULL)){

  p2 <- ggplot(res_dif, aes(x = age, y = diff)) + geom_line(size = 0.8) + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
          axis.line.x = element_line(colour = 'grey', size=0.2, linetype='solid'),
          axis.line.y = element_line(colour = 'grey', size=0.2, linetype='solid')) + 
    ylab("Difference men vs women")
  
  #plot with conf interval
  #ggplot(res_dif, aes(x = age, y = diff)) +
  #  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #  geom_line()
  }
  if (make_plots){
    print(p1)
    print(p2)
  }
  
  return (list("pdat" = pdat, "dif" = res_dif$diff, "inter_p" = gam.p,"g_beta" = gam.g_beta, "g_pv" = gam.g_pv))
} 




# Fit a GAM with age : gender interaction and (optional) correction for covariates (linear or spline)
fit_gam_age_only <- function(merged_tab, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots = T,  gam_family = gaussian(), min_age = 20, max_age = 80){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  
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
  return (pdat)
}


