library(RColorBrewer)
setwd("/Users/dashazhernakova/Documents/UMCG/data/gender_differences/omics/anova+gam")

# NMR metabolomics
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/Metabolomics/metabolomics_shared_folder/2.metabolites/LLD_bloodlipids_nmr.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- traits[,seq(5,ncol(traits))]

# proteins
traits <- as.data.frame(t(read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/rawProteinData/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
# drugs
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/LifeLines_phenotypes/20150618_45drugs_1135metasubjects_atleast5users.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
# pheno from Science
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/LifeLines_phenotypes/20150715_intristic_1135patients.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- traits[,seq(3,ncol(traits))]
# pheno more samples
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/LifeLines_phenotypes/20170123_selection_phenotypes_for_TL_not_binary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- traits[,seq(3,ncol(traits))]

cvd <- c("relativesCVD","hr","p","p_axis","pq","qrs","qrs_axis","qt","qtc","t_axis","dbp","hbf","map","sbp","cho","crph","glu","hb1c","hbac","hdc","k","ldc","tgl","ins","homair","bmi","hip","angioplasty","arythmia","chestpain","diabetes","diabtype","dilataorta","edema","heartattack","heartfailure","highchol","hypertension","lof","narrowcarotis","stroke","T1Dy_n","T2Dy_n","RelativesDiab","added_compl_legs_1_yes_0_no","added_pain_hands_or_feet_1y_0n","added_stiff_hands_or_feet_1y_0n","feetwounds","jointpainfeet","jointpainhands","legcomplnight","legcomplwalking","stiffnessfeet","stiffnesshands")
#traits <- as.data.frame(t(read.table(gzfile("/Users/dashazhernakova/Documents/UMCG/data/olink/LLDeep_expression/gene_level/LLD_genelevel_readcounts.nozeros.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

num_traits <- ncol(traits)

pheno <- as.data.frame(t(read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/correctedForAgeGenderSmokingContracCellCounts/age_gender_smk_contrac_cell_counts.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

match_idx = intersect(row.names(traits),row.names(pheno))
traits_m <- traits[match_idx,]
pheno_m <- pheno[match_idx,]
traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

#pheno_m$gender_F1M2 <- as.factor(pheno_m$gender_F1M2)
col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}


ntraits <- ncol(traits_m)
ntraits
nplotspp <- 20
npages = ceiling(ntraits/(2*nplotspp))
npages

cnt = 1
pdf("pheno_gam_dif.pdf", width = 15, height = 15)
par(mfrow=c(5,4))
for (idx in 1:ntraits){
  
  if (cnt > nplotspp){
    cnt = 1
    par(mfrow=c(5,4))
    print ("new page")
  }
  print(idx)
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  women <- merged_tab$gender_F1M2 == 1
  
  plot_scatter_and_gam2(merged_tab)
  cnt <- cnt + 1

}
dev.off()





rm_na_outliers <- function(traits_m, pheno_m, idx){
  traits_na.rm <- traits_m[!is.na(traits_m[,idx]),idx]
  pheno_na.rm <- pheno_m[!is.na(traits_m[,idx]),]
  women <- pheno_na.rm$gender_F1M2 == 1
  
  #merged_tab <- cbind(traits_na.rm, pheno_na.rm[,c("age", "gender_F1M2", "contrac")])
  #colnames(merged_tab) <- c(colnames(traits_m)[idx], "age", "gender_F1M2", "contrac")
  
  merged_tab <- cbind(traits_na.rm, pheno_na.rm)
  colnames(merged_tab) <- c(colnames(traits_m)[idx], colnames(pheno_na.rm))
  
  row.names(merged_tab) <- row.names(pheno_na.rm)
  
  # remove outliers further than 3 sds from the mean ( skip this if outcome is a factor )
  if (length(unique(merged_tab[,1])) > 3){
    w <- merged_tab[women,]
    m <- merged_tab[!women,]
    
    tab_nooutliers <- rbind(w[abs(w[,1] - mean(w[,1]))/sd(w[,1]) < 3,], m[abs(m[,1] - mean(m[,1]))/sd(m[,1]) < 3,])
  } else {
    tab_nooutliers <- merged_tab
  }
  return(tab_nooutliers)
}

plot_scatter_and_gam2 <- function(merged_tab){
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  pheno_name = colnames(merged_tab)[1]
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  women <- merged_tab$gender_F1M2 == 1
  
  a <-anova(lm(merged_tab[,1] ~ merged_tab$age + merged_tab$gender_F1M2 + merged_tab$age*merged_tab$gender_F1M2))
  anova_p <- a[3,"Pr(>F)"]
  
  if (length(unique(merged_tab[,1])) == 2){ # binary
    m1 <- gam(phenotype ~ gender_F1M2 + s(age, by = gender_F1M2), data = merged_tab, family=binomial(link = "logit"), method = "REML")
  } else {
    m1 <- gam(phenotype ~ gender_F1M2 + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
  }
  #m2 <- gam(phenotype ~  s(age, gender_F1M2, bs = "fs"), data = merged_tab, method = "REML")
  
  #summary(m1)
  m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), data = merged_tab, method = 'REML')
  m_o_p <- summary(m_o)$s.pv[2]
  
  pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = 500),gender_F1M2 = c('1', '2')))
  pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
  ylims <- with(merged_tab, range(phenotype))
  
  ## draw base plot
  palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
  plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
       main = paste0(pheno_name, "\nlinear interaction p = ", format(anova_p, digits = 3), "\nGAM interaction p = ", format(m_o_p, digits = 3)), 
       cex = 0.6, xlab = "age", ylab = pheno_name)
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
  }
  #plot(m_o, shade = TRUE,  scale = 0, seWithMean = TRUE, select = 2, main = "Trend in men as compared to women" , shift = coef(m0)[1])
}  
  
  
  #plot(m, shade = TRUE, pages = 1, scale = 0)
  
  pdat <- expand.grid(age = seq(20, 75, length = 500),
                      gender_F1M2 = c('1', '2'))
  xp <- predict(m, newdata = pdat, type = 'lpmatrix')
  
  ## which cols of xp relate to splines of interest?
  c1 <- grepl('gender_F1M21', colnames(xp))
  c2 <- grepl('gender_F1M22', colnames(xp))
  ## which rows of xp relate to sites of interest?
  r1 <- with(pdat, gender_F1M2 == '1')
  r2 <- with(pdat, gender_F1M2 == '2')
  X <- xp[r1, ] - xp[r2, ]
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(m)
  se <- sqrt(rowSums((X %*% vcov(m)) * X))
  crit <- qt(.975, df.residual(m))
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  res_dif <- data.frame(pair = "women-men",
             diff = dif,
             se = se,
             upper = upr,
             lower = lwr,
             age = seq(20, 75, length = 500))
  #ggplot(res_dif, aes(x = age, y = diff)) +
  #  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  #  geom_line()
    
 plot(res_dif$age, res_dif$diff, type = 'l')   

 
 
 
   
  #uri <- 'https://gist.githubusercontent.com/gavinsimpson/eb4ff24fa9924a588e6ee60dfae8746f/raw/geochimica-metals.csv'
  #metals <- read_csv(uri, skip = 1, col_types = c('ciccd'))
  #metals <- mutate(metals, SiteCode = factor(SiteCode))
  #m2 <- gam(Hg ~ SiteCode + s(Date, by = SiteCode), data = metals)
  #pdat2 <- expand.grid(Date = seq(1860, 2000, length = 400),
  #                    SiteCode = c('FION', 'CHNA', 'NODH'))
  #xp2 <- predict(m, newdata = pdat, type = 'lpmatrix')
  ### which cols of xp relate to splines of interest?
  #c12 <- grepl('CHNA', colnames(xp))
  #c22 <- grepl('FION', colnames(xp))
  ### which rows of xp relate to sites of interest?
  #r12 <- with(pdat, SiteCode == 'CHNA')
  #r22 <- with(pdat, SiteCode == 'FION')
  
  
  
  
  
  traits_na.rm <- merged_tab[,1]
  avg_w <- mean(traits_na.rm[women])
  avg_m <- mean(traits_na.rm[!women])
  step = 10
  diff_mw <- matrix(nrow = 7, ncol = 4)
  mean_norm <- matrix(nrow = 6, ncol = 3)
  bin_num <- 1
  for (bin_st in seq(20, 70, step)){
    bin_end <- bin_st + step - 1
    
    bin_w <- merged_tab$age >= bin_st & merged_tab$age < bin_end & merged_tab$gender_F1M2 == 1
    bin_m <- merged_tab$age >= bin_st & merged_tab$age < bin_end & merged_tab$gender_F1M2 == 2
    
    lines(c(bin_st, bin_end), c(mean(traits_na.rm[bin_w]), mean(traits_na.rm[bin_w])), lwd = 3, col = 1)
    lines(c(bin_st, bin_end), c(mean(traits_na.rm[bin_m]), mean(traits_na.rm[bin_m])), lwd = 3, col = 2)
    
    diff_mw[bin_num,1] <- bin_st + step/2
    diff_mw[bin_num,2] <- (mean(traits_na.rm[bin_m]) - mean(traits_na.rm[bin_w]))/mean(traits_na.rm)
    diff_mw[bin_num,3] <- (mean(traits_na.rm[bin_m]) - mean(traits_na.rm[bin_w]))/mean(traits_na.rm[bin_w])
    diff_mw[bin_num,4] <- (mean(traits_na.rm[bin_m]) - mean(traits_na.rm[bin_w]))
    #mean_norm[bin_num,1] <- bin_st + step/2
    #mean_norm[bin_num,2] <- mean(traits_na.rm[bin_w])/avg_w
    #mean_norm[bin_num,3] <- mean(traits_na.rm[bin_m])/avg_m
    
    bin_num <- bin_num + 1
  }
  plot(diff_mw[,1], diff_mw[,2], pch = 16, ylab = "(mean in men - mean in women)/mean in all", xlab = "age bin", main = colnames(traits_m)[idx])
  #plot(mean_norm[,1], mean_norm[,2], pch = 16, xlab = "age bin", main = colnames(traits_m)[idx], col = 1, ylim = c(0, max(mean_norm[,2],mean_norm[,3])))
  #points(mean_norm[,1], mean_norm[,3], pch = 16, col = 2)
  
