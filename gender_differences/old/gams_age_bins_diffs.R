setwd("/Users/dashazhernakova/Documents/UMCG/data/gender_differences/omics/test_norm_bin")
library("fANCOVA")
library(RColorBrewer)
library(mgcv)
#traits <- as.data.frame(t(read.table("/Users/dashazhernakova/Docuwoments/UMCG/data/olink/olink_shared_folder/Data/rawtraitseinData/CVD3_olinkNormal_1447_LLDsamples_t_traitsNames.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

#phenotypes
#traits <- read.table("../phenotypes_nonbinary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#diseases
#traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150715_diseases_1135patients.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#drugs
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_45drugs_1135metasubjects_atleast5users.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#diet 
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_diet_1135metasubjects_log_transform_imputed.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#pheno from Science
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_intrinsicfactors_1135metasubjects.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#smoking
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_smoke_1135metasubjects_noMissing.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

traits <- traits[,seq(3,ncol(traits))]

pheno <- as.data.frame(t(read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/correctedForAgeGenderSmokingContracCellCounts/age_gender_smk_contrac_cell_counts.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
match_idx = intersect(row.names(traits),row.names(pheno))
traits_m <- traits[match_idx,]
pheno_m <- pheno[match_idx,]
traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

nrow(pheno_m)
#palette(brewer.pal(n = 8, name = "Set1"))


col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

#span_val = 0.75

ntraits <- ncol(traits_m)
npages = 3
nplotspp <- ceiling(2*ntraits/npages)
pdf("pheno_science_gam_agebins_diff3.pdf", width = 15, height = 15)
for (page in 1:npages){
  par(mfrow=c(ceiling(nplotspp/5),5))
  
  for (np in 1:nplotspp){
    idx <- nplotspp*(page - 1) + np
    if (idx > ntraits){
      break
    }
    print(idx)
    traits_na.rm <- traits_m[!is.na(traits_m[,idx]),idx]
    pheno_na.rm <- pheno_m[!is.na(traits_m[,idx]),]
    women <- pheno_na.rm$gender_F1M2 == 1
    
    l1 <- lm(traits_na.rm~pheno_na.rm$age)
    l2 <- lm(traits_na.rm~pheno_na.rm$age + pheno_na.rm$gender_F1M2)
    l3 <- lm(traits_na.rm~pheno_na.rm$age + pheno_na.rm$gender_F1M2 + pheno_na.rm$age*pheno_na.rm$gender_F1M2)
    a <- anova(l1,l2, l3)
    anova_p <- a[3,"Pr(>F)"]
    
    palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
    plot(pheno_na.rm$age, traits_na.rm, col = pheno_na.rm$gender_F1M2,  pch = 16, main = colnames(traits_m)[idx], cex = 0.6, xlab = "age", ylab = colnames(traits_m)[idx])
    palette(c("indianred1", "dodgerblue1"))
    
    #gammod <- gam(traits_na.rm[women]~pheno_na.rm[women,"age"], family=binomial(link = "logit"))
    gammod <- gam(traits_na.rm[women]~s(pheno_na.rm[women,"age"]), method = "REML")
    lines(pheno_na.rm[women,"age"],fitted(gammod),col="indianred1", lwd = 2, lty = "dashed")
    gammod2 <- gam(traits_na.rm[!women]~s(pheno_na.rm[!women,"age"]), method = "REML")
    lines(pheno_na.rm[!women,"age"],fitted(gammod2),col="royalblue1", lwd = 2, lty = "dashed")
    
    avg_w <- mean(traits_na.rm[women])
    avg_m <- mean(traits_na.rm[!women])
    step = 10
    diff_mw <- matrix(nrow = 7, ncol = 2)
    mean_norm <- matrix(nrow = 6, ncol = 3)
    bin_num <- 1
    for (bin_st in seq(20, 70, step)){
      bin_end <- bin_st + step - 1
      
      bin_w <- pheno_na.rm$age >= bin_st & pheno_na.rm$age < bin_end & pheno_na.rm$gender_F1M2 == 1
      bin_m <- pheno_na.rm$age >= bin_st & pheno_na.rm$age < bin_end & pheno_na.rm$gender_F1M2 == 2
      
      lines(c(bin_st, bin_end), c(mean(traits_na.rm[bin_w]), mean(traits_na.rm[bin_w])), lwd = 3, col = 1)
      lines(c(bin_st, bin_end), c(mean(traits_na.rm[bin_m]), mean(traits_na.rm[bin_m])), lwd = 3, col = 2)
      
      diff_mw[bin_num,1] <- bin_st + step/2
      diff_mw[bin_num,2] <- (mean(traits_na.rm[bin_m]) - mean(traits_na.rm[bin_w]))/mean(traits_na.rm)
      
      #mean_norm[bin_num,1] <- bin_st + step/2
      #mean_norm[bin_num,2] <- mean(traits_na.rm[bin_w])/avg_w
      #mean_norm[bin_num,3] <- mean(traits_na.rm[bin_m])/avg_m
      
      bin_num <- bin_num + 1
    }
    plot(diff_mw[,1], diff_mw[,2], pch = 16, ylab = "(mean in men - mean in women)/mean in all", xlab = "age bin", main = colnames(traits_m)[idx], ylim = c(-0.5, 0.5))
    #plot(mean_norm[,1], mean_norm[,2], pch = 16, xlab = "age bin", main = colnames(traits_m)[idx], col = 1, ylim = c(0, max(mean_norm[,2],mean_norm[,3])))
    #points(mean_norm[,1], mean_norm[,3], pch = 16, col = 2)
    
  }
}
dev.off()

idx=5
idx=18
traits_na.rm <- traits_m[!is.na(traits_m[,idx]),idx]
pheno_na.rm <- pheno_m[!is.na(traits_m[,idx]),]
women <- pheno_na.rm$gender_F1M2 == 1

palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
plot(pheno_na.rm$age, traits_na.rm, col = pheno_na.rm$gender_F1M2,  pch = 16, main = colnames(traits_m)[idx], cex = 0.6, xlab = "age", ylab = colnames(traits_m)[idx])
palette(c("indianred1", "dodgerblue1"))

#gammod <- gam(traits_na.rm[women]~pheno_na.rm[women,"age"], family=binomial(link = "logit"))
gammod <- gam(traits_na.rm[women]~s(pheno_na.rm[women,"age"]), method = "REML")
lines(pheno_na.rm[women,"age"],fitted(gammod),col="indianred1", lwd = 2, lty = "dashed")
gammod2 <- gam(traits_na.rm[!women]~s(pheno_na.rm[!women,"age"]), method = "REML")
lines(pheno_na.rm[!women,"age"],fitted(gammod2),col="royalblue1", lwd = 2, lty = "dashed")

avg_w <- mean(traits_na.rm[women])
avg_m <- mean(traits_na.rm[!women])
step = 10
diff_mw <- matrix(nrow = 7, ncol = 2)
bin_num <- 1
for (bin_st in seq(20, 70, step)){
  bin_end <- bin_st + step - 1
  
  bin_w <- pheno_na.rm$age >= bin_st & pheno_na.rm$age < bin_end & pheno_na.rm$gender_F1M2 == 1
  bin_m <- pheno_na.rm$age >= bin_st & pheno_na.rm$age < bin_end & pheno_na.rm$gender_F1M2 == 2
  
  lines(c(bin_st, bin_end), c(mean(traits_na.rm[bin_w]), mean(traits_na.rm[bin_w])), lwd = 3, col = 1)
  lines(c(bin_st, bin_end), c(mean(traits_na.rm[bin_m]), mean(traits_na.rm[bin_m])), lwd = 3, col = 2)
  
  diff_mw[bin_num,1] <- bin_st + step/2
  diff_mw[bin_num,2] <- mean(traits_na.rm[bin_m]) - mean(traits_na.rm[bin_w])
  
  bin_num <- bin_num + 1
}
#plot(diff_mw[,1], diff_mw[,2], pch = 16, ylab = "mean in men - mean in women", xlab = "age bin", main = colnames(traits_m)[idx])


pheno_na.rm$age_bin <- trunc(pheno_na.rm$age/10)

res_aov <- aov(traits_na.rm~pheno_na.rm$age + pheno_na.rm$gender_F1M2 + pheno_na.rm$age*pheno_na.rm$gender_F1M2)
res_aov2 <- aov(traits_na.rm~pheno_na.rm$age_bin + pheno_na.rm$gender_F1M2 + pheno_na.rm$age_bin*pheno_na.rm$gender_F1M2)
summary(res_aov)
summary(res_aov2)
