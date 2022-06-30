args <- commandArgs(trailingOnly = TRUE)
library(RColorBrewer)
library(ggplot2)

source("/groups/umcg-lifelines/tmp01/projects/ov20_0051//umcg-dzhernakova/scripts/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
setwd("/groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/gender_difs/factors")


col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

plot_3d <- function(pheno_name, gam.fit, merged_tab, covariates, predictor,  out_path = "./", min_age = 20, max_age = 80, plot_points = T, n_points = 300){
  
  if (is.factor(merged_tab[,predictor])){
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = 300),
                                          gender_F1M2 = c('1', '2'),
                                          predictor = levels(merged_tab[,predictor])))
  } else {
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = 300),
                                          gender_F1M2 = c('1', '2'),
                                          predictor = seq(quantile(merged_tab[,predictor], probs = 0.01), quantile(merged_tab[,predictor], probs = 0.99), length = 50)))
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
  
  log_tr <- c("LEU","LY","MO","MOP","GR","BA","BAP","EO","EOP","TGL","HAL1","HALB","AST","ALT","AF","GGT","LCRP","TSH","UKRO","UKR24")
  if (pheno_name %in% log_tr){
    merged_tab[,1] <-  exp(merged_tab[,1]) - 1
    pdat$pred <- exp(pdat$pred) - 1
    pdat$lwr <- exp(pdat$lwr) - 1
    pdat$upr <- exp(pdat$upr) - 1
  }
  
  
  #pheno_name = "CHO"
  colnames(pdat) <- gsub(predictor,"predictor", colnames(pdat))
  pdf(paste0(out_path, "/", pheno_name, "_vs_", predictor, ".3d.pdf"), height = 5, width = 10)
  p <- ggplot(pdat, aes(x = age, y = predictor,fill=pred )) +
    geom_raster() +
    geom_contour(aes_(z = ~ pred), bins = 10, colour = "grey") +
    facet_wrap(~gender_F1M2, ncol =2) +
    geom_hline(yintercept = quantile(merged_tab[,predictor], probs = 0.05), color = "orange1") +
    geom_hline(yintercept = quantile(merged_tab[,predictor], probs = 0.95), color = "cadetblue1") +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    labs(x="age", y = predictor, fill = pheno_name, title = paste0(pheno_name, " by age and ",predictor, ", sex:")) +         theme_minimal()
  print(p)
  #pdf("test_3d.pdf", height = 5, width =10)
  #dev.off()
  #return (p)
  dev.off()
}


plot_2d <- function(pheno_name, gam.fit, merged_tab, covariates, predictor,  out_path = "./", min_age = 20, max_age = 80, plot_points = T, n_points = 300){
  
  if (is.factor(merged_tab[,predictor])){
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = 300),
                                          gender_F1M2 = c('1', '2'),
                                          predictor = levels(merged_tab[,predictor])))
  } else {
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = 300),
                                          gender_F1M2 = c('1', '2'),
                                          predictor = seq(quantile(merged_tab[,predictor], probs = 0.05), quantile(merged_tab[,predictor], probs = 0.95), length = 5)))
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
  
  log_tr <- c("LEU","LY","MO","MOP","GR","BA","BAP","EO","EOP","TGL","HAL1","HALB","AST","ALT","AF","GGT","LCRP","TSH","UKRO","UKR24")
  if (pheno_name %in% log_tr){
    merged_tab[,1] <-  exp(merged_tab[,1]) - 1
    pdat$pred <- exp(pdat$pred) - 1
    pdat$lwr <- exp(pdat$lwr) - 1
    pdat$upr <- exp(pdat$upr) - 1
  }
  
  
  #pheno_name = "CHO"
  colnames(pdat) <- gsub(predictor,"predictor", colnames(pdat))
  pdf(paste0(out_path, "/", pheno_name, "_vs_", predictor, ".2d.pdf"), height = 5, width = 10)
  p <- ggplot(pdat, aes(x = age, y = pred, color= as.factor(predictor) )) +
    facet_wrap(~gender_F1M2, ncol =2) +
    geom_line() +
    scale_color_brewer(palette = "Spectral", direction = -1) +
    labs(x="age", y = pheno_name,  color = predictor, title = paste0(pheno_name, " by age and ",predictor, ", sex:")) +
    theme_minimal()
  print(p)
  #pdf("test_3d.pdf", height = 5, width =10)
  #dev.off()
  #return (p)
  dev.off()
}

plot_res <- function(pheno_name, gam.fit, merged_tab, covariates, predictor, out_path = "./", min_age = 20, max_age = 80, plot_points = F, n_points = 300){
  
  if (is.factor(merged_tab[,predictor])){
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = 300),
                                          gender_F1M2 = c('1', '2'),
                                          predictor = levels(merged_tab[,predictor])))
  } else {
    new.x <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = 300),
                                          gender_F1M2 = c('1', '2'),
                                          predictor = c(quantile(merged_tab[,predictor], probs = 0.05), quantile(merged_tab[,predictor], probs = 0.95))))
  }
  cat(pheno_name, quantile(merged_tab[,predictor], probs = 0.05), quantile(merged_tab[,predictor], probs = 0.95), "\n" )
  
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
  log_tr <- c("LEU","LY","MO","MOP","GR","BA","BAP","EO","EOP","TGL","HAL1","HALB","AST","ALT","AF","GGT","LCRP","TSH","UKRO","UKR24")
  if (pheno_name %in% log_tr){
    merged_tab[,1] <-  exp(merged_tab[,1]) - 1
    pdat$pred <- exp(pdat$pred) - 1
    pdat$lwr <- exp(pdat$lwr) - 1
    pdat$upr <- exp(pdat$upr) - 1
  }
  #pheno_name = "CHO"
  cex_main = 1
  plot_title = paste0(pheno_name, " plotted by", predictor)
  
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
  pdf(paste0(out_path, "/", pheno_name, "_plotby_", predictor, ".pdf"))
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

plot_res_by_factor <- function(pheno_name, gam.fit, merged_tab, covariates, predictor, min_age = -30, max_age = 30, plot_points = T, n_points = 300){
  
  #if (is.factor(merged_tab[,predictor])){
  #  new.x <- with(merged_tab, expand.grid(age = seq(-30, 30, length = 300),
  #                                        gender_F1M2 = c('1', '2'),
  #                                        predictor = levels(merged_tab[,predictor])))
  #} else {
  min_val <- min(merged_tab[,predictor])
  max_val <- max(merged_tab[,predictor])
  new.x <- with(merged_tab, expand.grid(age = c(-27, 27),
                                        gender_F1M2 = c('1', '2'),
                                        predictor = seq(min_val, max_val, length = 300)))
  #}
  
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
  
  pdat[,"age"] <- as.factor(pdat[,"age"])
  levels( pdat[,"age"]) <- c("low", "high")
  
  #pheno_name = "CHO"
  cex_main = 1
  plot_title = paste0(pheno_name, " plotted by", predictor)
  
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
  pdf(paste0(out_path, "/", pheno_name, "_vs_", predictor, "_byage.pdf"))
  palette(c(col2transparent("#ff9999", 120),col2transparent("#99ccff", 120)))
  par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i") # Remove plot padding
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  
  if (plot_points){
    plot(merged_tab2$phenotype ~ merged_tab2[, predictor], col = merged_tab2$gender_F1M2,  pch = 16,
         main = plot_title,
         cex = 0.6, xlab = predictor, ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T,
         ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
         xlim = c(min_val, max_val))
    
    
  } else{
    plot(1, type="n",
         main = plot_title,
         cex = 0.6, xlab = predictor, ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T,
         ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
         xlim = c(min(merged_tab2$age), max(merged_tab2$age)))
  }
  
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2[, predictor]), col = "grey90")
  
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
    dd_low <- dd[dd[,"age"] == "low", ]
    lines(dd_low$pred ~ dd_low[,predictor],  col = cols[[l]], lwd = 2, lty = 1)
    polygon(c(rev(dd_low[,predictor]), dd_low[,predictor]), c(rev(dd_low$lwr), dd_low$upr), col = col2transparent(cols[[l]], 65), border = NA)
    
    dd_high <- dd[dd[,"age"] == "high", ]
    lines(dd_high$pred ~ dd_high[,predictor],  col = cols2[[l]], lwd = 2, lty = 1)
    polygon(c(rev(dd_high[,predictor]), dd_high[,predictor]), c(rev(dd_high$lwr), dd_high$upr), col = col2transparent(cols[[l]], 65), border = NA)
    
    
  }
  dev.off()
  
}

rm_outliers <- function(merged_tab){
  merged_tab <- na.omit(merged_tab)
  merged_tab <- merged_tab[(merged_tab$age < 80) & (merged_tab$age >= 20),]
  m <- merged_tab
  mq1 <- quantile(m[,1], probs = 0.25)
  mq3 <- quantile(m[,1], probs = 0.75)
  miqr <- mq3 - mq1
  m_clean <- m[m[,1] < mq3 + 3*miqr & m[,1] > mq1 - 3*miqr,]
  
  return(m_clean)
}

traits_path <- "../v5/data/LL_phenotypes_merged_all.log_some.v5.txt"
pheno_path <- "factors_final.v2.txt"
pheno_name=args[1]
med=args[2]
if (length(args) > 2){
  run_bootstrap=args[3]
} else {
  run_bootstrap=F
}
#pheno_name="CHO"
#med="statins"

out_path<-"results_v5/"
cat("Data paths:\nphenotype traits:", traits_path, "\r\ncovariates:", pheno_path, "phenotype:", pheno_name, "medication name:", med, "\n")
print("v3!")
# read phenotype traits of interest
traits <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

# read age, gender and other covariate phenotypes
pheno <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)


#order samples in the two tables
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)

cols2use <- c("gender_F1M2", "age", "smoking", "stress_chronic", "phys_activity_total", "diet", "alcohol", "sleep_duration")
if (med != "NA"){
  cat("Medication will be used as a covariate\n")
  cols2use <- c(cols2use, med)
  merged_tab <- cbind(traits_m[, pheno_name], pheno_m[,cols2use])
  colnames(merged_tab)[1] <- "phenotype"
  colnames(merged_tab)[ncol(merged_tab)] <- "med"
  
  merged_tab <- rm_outliers(merged_tab)
  
  # preds <-  colnames(merged_tab)[c(4,5,8:(ncol(merged_tab)-1))]
  preds <-   colnames(merged_tab)[5:(ncol(merged_tab) -1)]
  terms_binary <- " + smoking + s(age, by = smoking) + s(age, by = interaction(smoking, gender_F1M2)) + med + s(age, by = med) + s(age, by = interaction(med, gender_F1M2))"
  merged_tab$med <-as.factor(merged_tab$med)
  
} else {
  cat("No medication to use\n")
  merged_tab <- cbind(traits_m[, pheno_name], pheno_m[,cols2use])
  colnames(merged_tab)[1] <- "phenotype"
  
  merged_tab <- rm_outliers(merged_tab)
  
  preds <-  colnames(merged_tab)[5:ncol(merged_tab)]
  terms_binary <- " + smoking + s(age, by = smoking) + s(age, by = interaction(smoking, gender_F1M2))"
  
}
#center age at mean age
#merged_tab$age <- merged_tab$age - 50

merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
merged_tab$smoking <- gsub("2", "0", merged_tab$smoking)
merged_tab$smoking <- as.factor(merged_tab$smoking)
#merged_tab <- mutate(merged_tab, smoking = ordered(smoking, levels = c('0', '1')))


terms <-  paste0(" + s(", paste(preds, collapse = ")+ s("), ")")
terms_inter3 <- paste0(" + ti(", paste(preds, collapse = ", age, by = gender_F1M2)+ ti("), ", age, by = gender_F1M2)")
terms_inter_age <- paste0(" + ti(", paste(preds, collapse = ", age)+ ti("), ", age)")
terms_inter_sex <- paste0(" + s(", paste(preds, collapse = ", by = gender_F1M2)+ s("), ", by = gender_F1M2)")


full_formula <- as.formula(paste("phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2) ", terms, terms_inter3, terms_inter_age, terms_inter_sex, terms_binary, sep = " "))
print(full_formula)

if (run_bootstrap == "T"){
  cat("Running bootstrap!\n")
  ids <- sample(nrow(merged_tab), nrow(merged_tab), replace = T)
  merged_tab <- merged_tab[ids,]
}


full_fit <- gam(full_formula, data = merged_tab, method = "REML", select=T)

#fit <- gam(phenotype ~ smoking + SMK3 + gender_F1M2 + s(age) + s(LTE_SUM) + s(LDI_SUM) + s(total_mwk_VAL) + s(  total_scor_VAL) + s(MVPA_mwk_VAL) + s(MVPA_scor_VAL) +  s(LLDS_T1A) + s(SumOfalcohol),data = merged_tab, method = "REML", select=T)
s <- summary(full_fit)
print(s)

if (run_bootstrap == "T"){
  pheno_name=paste0(pheno_name, "_bootstrap_", args[4])
}
write.table(s$p.table, file = paste0(out_path, pheno_name, ".p.table.txt"), sep = "\t", quote = F, col.names = NA)
write.table(s$s.table, file = paste0(out_path, pheno_name, ".s.table.txt"), sep = "\t", quote = F, col.names = NA)
#}



#
# Plotting
#
if (run_bootstrap != "T"){
  if (med != "NA"){
    covariates <- c(preds, "smoking", "med")
  } else {
    covariates <- c(preds, "smoking")
  }
  
  for (predictor in covariates){
    plot_3d(pheno_name, full_fit, merged_tab, covariates, predictor,out_path)
    plot_2d(pheno_name, full_fit, merged_tab, covariates, predictor,out_path)
    plot_res(pheno_name, full_fit, merged_tab, covariates, predictor, out_path)
  }
}
