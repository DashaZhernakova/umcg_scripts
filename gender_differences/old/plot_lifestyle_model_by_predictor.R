library(RColorBrewer)

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}
plot_res <- function(pheno_name, gam.fit, merged_tab, covariates, predictor, min_age = -30, max_age = 30, plot_points = T, n_points = 300){
  
  if (is.factor(merged_tab[,predictor])){
    new.x <- with(merged_tab, expand.grid(age = seq(-30, 30, length = 300), 
                                          gender_F1M2 = c('1', '2'),
                                          predictor = levels(merged_tab[,predictor]))) 
  } else {
    new.x <- with(merged_tab, expand.grid(age = seq(-30, 30, length = 300), 
                                          gender_F1M2 = c('1', '2'),
                                          predictor = c(quantile(merged_tab[,predictor], probs = 0.05), quantile(merged_tab[,predictor], probs = 0.95))))
  }
  
  for (c in covariates){
    if (c != predictor){
      if (is.factor(merged_tab[,c])){
        new.x[,c] <- 0
        new.x[,c] <- as.factor(new.x[,c])
      } else {
        new.x[,c] <- mean(merged_tab[,c])
      }
    }
  }
  colnames(new.x) <- gsub("predictor", predictor,colnames(new.x))
  new.y <- data.frame(predict(gam.fit, newdata = new.x, se.fit = TRUE, type = "response"))
  pdat <- data.frame(new.x, new.y)
  pdat <- rename(pdat, pred = fit, SE = se.fit)
  pdat <- mutate(pdat, lwr = pred - 1.96 * SE, upr = pred + 1.96 * SE) # calculating the 95% confidence interval
  
  pdat[,predictor] <- as.factor(pdat[,predictor])
  levels( pdat[,predictor]) <- c("low", "high")
  
  #pheno_name = "CHO"
  cex_main = 1
  plot_title = pheno_name
  
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
  #palette(c(col2transparent("indianred1", alpha_points),col2transparent("dodgerblue1", alpha_points)))
  pdf(paste0(pheno_name, "_plotby_", predictor, ".pdf"))
  palette(c(col2transparent("#ff9999", 120),col2transparent("#99ccff", 120)))
  par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i") # Remove plot padding
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  
    if (plot_points){
      plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
           main = plot_title, 
           cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
           ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
           xlim = c(min(merged_tab2$age), max(merged_tab2$age)))
      
      #if (length(factor_name) > 0){
      #  points(merged_tab2[merged_tab2[,factor_name] == 1,"age"], merged_tab2[merged_tab2[,factor_name],"phenotype"], pch = 8, col = col2transparent("gold", 120), cex = 0.6)
      #}
    } else{
      plot(1, type="n", 
           main = plot_title, 
           cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
           ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
           xlim = c(min(merged_tab2$age), max(merged_tab2$age)))
    }
  
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2$age), col = "grey90")
  
  #at = pretty(merged_tab2$age)
  #mtext(side = 1, text = at, at = at, 
  #      col = "grey20", line = 1, cex = 0.4)
  
  #at = pretty(merged_tab2$phenotype)  
  #mtext(side = 2, text = at, at = at, col = "grey20", line = 1, cex = 0.4)
  
   levs <- levels(merged_tab$gender_F1M2)
  cols = list("indianred1", "dodgerblue1")
  cols2 <- list("orange1","cadetblue1")
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    dd_low <- dd[dd[,predictor] == "low", ]
    lines(pred ~ age, data = dd_low, col = cols[[l]], lwd = 2, lty = 1)
    polygon(c(rev(dd_low$age), dd_low$age), c(rev(dd_low$lwr), dd_low$upr), col = col2transparent(cols[[l]], 65), border = NA)
    
    dd_high <- dd[dd[,predictor] == "high", ]
    lines(pred ~ age, data = dd_high, col = cols2[[l]], lwd = 2, lty = 1)
    polygon(c(rev(dd_high$age), dd_high$age), c(rev(dd_high$lwr), dd_high$upr), col = col2transparent(cols[[l]], 65), border = NA)
    
    
  }
  dev.off()
  
}

# source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
# setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/factors")
# 
# rm_outliers <- function(merged_tab){
#   merged_tab <- na.omit(merged_tab)
#   merged_tab <- merged_tab[(merged_tab$age < 80) & (merged_tab$age >= 20),]
#   
#   w <- merged_tab[merged_tab$gender_F1M2 == 1,]
#   m <- merged_tab[merged_tab$gender_F1M2 == 2,]
#   mq1 <- quantile(m[,1], probs = 0.25)
#   mq3 <- quantile(m[,1], probs = 0.75)
#   miqr <- mq3 - mq1
#   m_clean <- m[m[,1] < mq3 + 1.5*miqr & m[,1] > mq1 - 1.5*miqr,]
#   
#   wq1 <- quantile(w[,1], probs = 0.25)
#   wq3 <- quantile(w[,1], probs = 0.75)
#   wiqr <- wq3 - wq1
#   w_clean <- w[w[,1] < wq3 + 1.5*wiqr & w[,1] > wq1 - 1.5*wiqr,]
#   
#   tab_nooutliers <- rbind(w_clean, m_clean)
#   return(tab_nooutliers)
# }
# 
# 
# traits_path <- "../v4/data/LL_phenotypes_merged_all.log_some.v5.txt"
# pheno_path <- "factors_final.txt"
# pheno_name="CHO"
# med="statins"
# out_path<-"/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/factors/results/"
# cat("Data paths:\nphenotype traits:", traits_path, "\r\ncovariates:", pheno_path, "phenotype:", pheno_name, "medication name:", med, "\n")
# print("v2!")
# # read phenotype traits of interest
# traits <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
# 
# # read age, gender and other covariate phenotypes
# pheno <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
# 
# 
# #order samples in the two tables
# traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
# pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
# all(row.names(traits_m) == row.names(pheno_m))
# num_traits <- ncol(traits_m)
# merged_tab <- cbind(traits_m[, pheno_name], pheno_m[,1:9], pheno_m[,med])
# colnames(merged_tab)[1] <- "phenotype"
# colnames(merged_tab)[ncol(merged_tab)] <- "med"
# 
# merged_tab <- rm_outliers(merged_tab)
# #center age at mean age
# mean_age = 50
# merged_tab$age <- merged_tab$age - mean_age
# 
# # preds <-  colnames(merged_tab)[c(4,5,8:(ncol(merged_tab)-1))]
# preds <-  colnames(merged_tab)[c(4,5,7:(ncol(merged_tab)-1))]
# 
# merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
# merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
# m <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2) + s(phys_activity_intensive) + ti(phys_activity_intensive, age, by = ord_gender_F1M2) + s(phys_activity_intensive, by = ord_gender_F1M2), data = merged_tab, method = "REML", select=T)
# gam.fit <- m
# covariates <- character(0)
# predictor <- "phys_activity_intensive"



