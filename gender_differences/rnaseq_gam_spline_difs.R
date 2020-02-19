library(RColorBrewer)
library('dplyr')
library('mgcv')
#setwd("/Users/dashazhernakova/Documents/UMCG/data/gender_differences/omics/anova+gam")
setwd('/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/')



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
  pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  women <- merged_tab$gender_F1M2 == 1
  
  lm1 <- lm(merged_tab[,1] ~ merged_tab$age + merged_tab$gender_F1M2 + merged_tab$age*merged_tab$gender_F1M2)
  #a <-anova(lm1)
  #anova_p <- a[3,"Pr(>F)"]
  
  if (length(unique(merged_tab[,1])) == 2){ # binary
    m1 <- gam(phenotype ~ gender_F1M2 + s(age, by = gender_F1M2), data = merged_tab, family=binomial(link = "logit"), method = "REML")
  } else {
    m1 <- gam(phenotype ~ gender_F1M2 + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
  }
  #m2 <- gam(phenotype ~  s(age, gender_F1M2, bs = "fs"), data = merged_tab, method = "REML")
  
  #summary(m1)
  m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), data = merged_tab, method = 'REML')
  m_o_p <- summary(m_o)$s.pv[2]
  if (m_o_p > 0.05){
    return (NULL)
  }
  pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = 500),gender_F1M2 = c('1', '2')))
  pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
  ylims <- with(merged_tab, range(phenotype))
  
  ## draw base plot
  palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
  plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
       main = paste0(pheno_name, "\nGAM interaction p = ", format(m_o_p, digits = 3)), 
       cex = 0.6, xlab = "age", ylab = pheno_name)
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
  }
  #plot(m_o, shade = TRUE,  scale = 0, seWithMean = TRUE, select = 2, main = "Trend in men as compared to women" , shift = coef(m_o)[1])
  
  pdat <- expand.grid(age = seq(20, 75, length = 400),
                      gender_F1M2 = c('1', '2'))
  res_dif <- smooth_diff(m1, pdat, '2', '1', "gender_F1M2")
  res_dif$age = seq(20, 75, length = 400)
  plot(res_dif$age, res_dif$diff, type = 'l')
  return (res_dif)
}  



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



# expression
traits0 <- as.data.frame(t(read.table(gzfile("/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/LLD_genelevel_readcounts.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
gte <- read.table("/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/gte_all.txt", sep = "\t", as.is = T, check.names = F)
gte_m <- gte[match(row.names(traits0), gte[,2], nomatch = 0 ),]
traits <- traits0[match(gte_m[,2], row.names(traits0), nomatch = 0 ),]
all(row.names(traits) == gte_m[,2])
row.names(traits) = gte_m[,1]
#cvd <- c("relativesCVD","hr","p","p_axis","pq","qrs","qrs_axis","qt","qtc","t_axis","dbp","hbf","map","sbp","cho","crph","glu","hb1c","hbac","hdc","k","ldc","tgl","ins","homair","bmi","hip","angioplasty","arythmia","chestpain","diabetes","diabtype","dilataorta","edema","heartattack","heartfailure","highchol","hypertension","lof","narrowcarotis","stroke","T1Dy_n","T2Dy_n","RelativesDiab","added_compl_legs_1_yes_0_no","added_pain_hands_or_feet_1y_0n","added_stiff_hands_or_feet_1y_0n","feetwounds","jointpainfeet","jointpainhands","legcomplnight","legcomplwalking","stiffnessfeet","stiffnesshands")
#traits <- as.data.frame(t(read.table(gzfile("/Users/dashazhernakova/Documents/UMCG/data/olink/LLDeep_expression/gene_level/LLD_genelevel_readcounts.nozeros.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

num_traits <- ncol(traits)

pheno <- as.data.frame(t(read.table("/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/age_gender_smk_contrac_cell_counts.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]



gene_table <- read.table("/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/geneid_to_gene.txt", header = T, sep = "\t", as.is = T, check.names = F)


#pheno_m$gender_F1M2 <- as.factor(pheno_m$gender_F1M2)
col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}


ntraits <- ncol(traits_m)
ntraits
nplotspp <-20
npages = ceiling(2*ntraits/(nplotspp))
npages

res_dif_all <- data.frame(age = seq(20, 75, length = 400))

cnt = 1
pdf("expression_all_gam_dif_sign.pdf", width = 15, height = 15)
par(mfrow=c(5,4))
for (idx in 1:ntraits){
  
  if (cnt > nplotspp){
    cnt = 1
    par(mfrow=c(5,4))
    print ("new page")
  }
  print(idx)
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  
  res_dif = NULL
  res_dif <- plot_scatter_and_gam2(merged_tab)
  if (! is.null(res_dif)){
    break
    cnt <- cnt + 1
    #  pheno <- colnames(merged_tab)[1]
    #  res_dif_all[,pheno] <- res_dif$diff
  }
  
  
}
dev.off()
