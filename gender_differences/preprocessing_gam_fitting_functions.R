library(RColorBrewer)
library('dplyr')
library('mgcv')


# Combine phenotype and trait of interest and remove outliers in the trai of interest
rm_na_outliers <- function(traits_m, pheno_m, idx, method = "IQR", scale_tr = F, log_tr = F, int_tr = F){
  traits_na.rm <- traits_m[!is.na(traits_m[,idx]),idx]
  pheno_na.rm <- pheno_m[!is.na(traits_m[,idx]),]
  women <- pheno_na.rm$gender_F1M2 == 1
  
  #merged_tab <- cbind(traits_na.rm, pheno_na.rm[,c("age", "gender_F1M2", "contrac")])
  #colnames(merged_tab) <- c(colnames(traits_m)[idx], "age", "gender_F1M2", "contrac")
  all(row.names(traits_na.rm) == row.names(pheno_na.rm))
  
  merged_tab <- cbind(traits_na.rm, pheno_na.rm)
  colnames(merged_tab) <- c(colnames(traits_m)[idx], colnames(pheno_na.rm))
  
  row.names(merged_tab) <- row.names(pheno_na.rm)
  
  # remove outliers ( skip this if outcome is a factor )
  if (length(unique(merged_tab[,1])) > 3){
    w <- merged_tab[women,]
    m <- merged_tab[!women,]
   
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
      print ("Wrong method! No outlier removal!")
      tab_nooutliers <- merged_tab
    }
    
  } else {
    tab_nooutliers <- merged_tab
  }
  
  if (log_tr & scale_tr){ 
    tab_nooutliers[,1] <- scale(log(tab_nooutliers[,1] + 1))
  } else if (log_tr & ! scale_tr) {
    tab_nooutliers[,1] <- log(tab_nooutliers[,1] + 1)
  } else if (!log_tr & scale_tr) {
    tab_nooutliers[,1] <- scale(tab_nooutliers[,1])
  } else if (int_tr) {
    tab_nooutliers[,1] <- qnorm((rank(tab_nooutliers[,1],na.last="keep")-0.5)/sum(!is.na(tab_nooutliers[,1])))
  }
  return(tab_nooutliers)
}

# Fit a GAM with age : gender interaction and (optional) correction for cell counts
plot_scatter_and_gam2 <- function(merged_tab, pheno_name, correctForCellCounts, n_points = 300, make_plots, label = '', gam_family = gaussian()){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  #pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  
  women <- merged_tab$gender_F1M2 == 1
  res_dif = NULL
  
  if (!correctForCellCounts){
    
      m1 <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML", family = gam_family)
      m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
                 data = merged_tab, method = 'REML', family = gam_family)
    
    
    m_o_p <- summary(m_o)$s.pv[length(summary(m_o)$s.pv)]
    m_o_g_beta <- m_o$coefficients["ord_gender_F1M2.L"]
    m_o_g_pv <- summary(m_o)$p.pv["ord_gender_F1M2.L"]
    
    if (m_o_p > 0.05){
      return (list("dif" = NULL, "inter_p" = m_o_p,"g_beta" = m_o_g_beta, "g_pv" = m_o_g_pv))
    }
    pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                         gender_F1M2 = c('1', '2')))
    pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
    
    
    
    if (! "Ordered Categorical" %in% gam_family$family){
      pdat_dif <- expand.grid(age = seq(20, 75, length = n_points),
                              gender_F1M2 = c('1', '2'))
      #res_dif <- smooth_diff(m1, pdat_dif, '2', '1', "gender_F1M2")
      #res_dif$age = seq(20, 75, length = n_points)
      res_dif <- simple_diff(pdat)
    } else { # stupid way to go
      m1 <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
      pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                           gender_F1M2 = c('1', '2')))
      pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
      
    }
  } else { # Correct for cell counts
    
    m1 <- gam(phenotype ~ gender_F1M2 + ba + eo + er + gr + 
                ly + mo + tr + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
    
    m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + ba + eo + er + gr + 
                 ly + mo + tr + s(age, by = ord_gender_F1M2), 
               data = merged_tab, method = 'REML')
    m_o_p <- summary(m_o)$s.pv[length(summary(m_o)$s.pv)]
    m_o_g_beta <- m_o$coefficients["ord_gender_F1M2.L"]
    m_o_g_pv <- summary(m_o)$p.pv["ord_gender_F1M2.L"]
    
    if (m_o_p > 0.05){
      return (list("dif" = NULL, "inter_p" = m_o_p,"g_beta" = m_o_g_beta, "g_pv" = m_o_g_pv))
    }
    
    pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), gender_F1M2 = c('1', '2'), 
                                         ba = mean(ba), eo = mean(eo), er = mean(er), gr = mean(gr), 
                                         ly = mean(ly),  mo = mean(mo), tr = mean(tr)))
    pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
    
    if (! "Ordered Categorical" %in% gam_family$family){
      pdat_dif <- expand.grid(age = seq(20, 75, length = n_points), gender_F1M2 = c('1', '2'), 
                              ba = mean(merged_tab$ba), eo = mean(merged_tab$eo), er = mean(merged_tab$er), gr = mean(merged_tab$gr), 
                              ly = mean(merged_tab$ly),  mo = mean(merged_tab$mo), tr = mean(merged_tab$tr))
      #res_dif <- smooth_diff(m1, pdat_dif, '2', '1', "gender_F1M2")
      #res_dif$age = seq(20, 75, length = n_points)
      res_dif <- simple_diff(pdat)
    } else { # stupid way to go
      m1 <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
      pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                           gender_F1M2 = c('1', '2')))
      pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
    }
  }
  
  
  ylims <- with(merged_tab, range(phenotype))
  ## draw base plot
  palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
  #p1 <- plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
  #     main = paste0(pheno_name, ' ', label, "\nGAM interaction p = ", format(m_o_p, digits = 3)), 
  #     cex = 0.6, xlab = "age", ylab = pheno_name)
  
  p1 <- ggplot(merged_tab, aes(y = phenotype, x = age,  col = gender_F1M2)) + geom_point(aes(col = gender_F1M2)) + 
    scale_color_manual(values = c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125),"indianred1", "dodgerblue1")) + 
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5), 
          axis.line.x = element_line(colour = 'grey', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'grey', size=0.5, linetype='solid')) + 
    ylab(pheno_name) + 
    ggtitle(paste0(pheno_name, ' ', label, "\nGAM interaction p = ", format(m_o_p, digits = 3))) + 
    scale_x_discrete(name ="age", limits=seq(10,80,10))
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    #lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
    p1 <- p1 + geom_line(data = dd, aes(x= age, y = pred), colour = cols[[l]], size = 1)
  }
  
  #Plot the difference
  if (! is.null(res_dif = NULL)){
  #p2 <- plot(res_dif$age, res_dif$diff, type = 'l')
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
    cex_main = 1
    ylims <- with(merged_tab, range(phenotype))
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
         main = paste0(pheno_name, ' ', label, "\nGAM interaction p = ", format(m_o_p, digits = 3)), 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main)
    
    levs <- levels(merged_tab$gender_F1M2)
    cols = c("indianred1", "dodgerblue1")
    
    ## add the fitted lines
    for (l in seq_along(levs)) {
      dd <- pdat[pdat$gender_F1M2 == levs[l],]
      lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
      
    }
    if (! is.null(res_dif = NULL)){
      #Plot the difference
      plot(res_dif$age, res_dif$diff, type = 'l')
    }
  }
  
  return (list("pdat" = pdat, "dif" = res_dif$diff, "inter_p" = m_o_p,"g_beta" = m_o_g_beta, "g_pv" = m_o_g_pv))
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

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}


test_polynomial_interaction <- function(merged_tab, model_fam = NO(), gm_mean_cov = F){
  if (! gm_mean_cov){
    #lm1 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) + I(age^4), data = merged_tab, family = model_fam, method = RS(30))
    #lm2 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^3) + I(age^4) + gender_F1M2:(age + I(age^2) + I(age^3) + I(age^4)), data = merged_tab, family = model_fam, method = RS(30)) 
    lm1 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) , data = merged_tab, family = model_fam, method = RS(30))
    lm2 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)), data = merged_tab, family = model_fam, method = RS(30)) 
  } else {
    lm1 <- gamlss(phenotype ~ gender_F1M2 + age + gm_mean + I(age^2) + I(age^3) , data = merged_tab, family = model_fam, method = RS(30))
    lm2 <- gamlss(phenotype ~ gender_F1M2 + age + gm_mean + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)), data = merged_tab, family = model_fam, method = RS(30)) 
    
  }
  
  lr <- LR.test(lm1, lm2, print = F)
  m <- lm1
  if (lr$p.val < 0.05){
    m <- lm2
  }
  return (list("lr_p.val" = lr$p.val, "mod" = m))
}

plot_scatter_and_gamlss <- function(merged_tab, pheno_name, n_points = 300, make_plots = T, gam_family = NO(), label = ''){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  
  res <- test_polynomial_interaction(merged_tab, gam_family)  
  inter_pval <- res[["lr_p.val"]]
  m1 <- res[["mod"]] 
  if (inter_pval > 0.05){
      return (NULL)
  }

  pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                       gender_F1M2 = c('1', '2')))
  pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
  pdat_dif <- expand.grid(age = seq(20, 75, length = n_points),
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
