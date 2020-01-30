library("fANCOVA")
library(RColorBrewer)
library(mgcv)

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

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}


ntraits <- ncol(traits_m)
ntraits
nplotspp <- 20
npages = ceiling(ntraits/nplotspp)
npages

cnt = 1
pdf("drugs_anova_gam.pdf", width = 15, height = 15)
par(mfrow=c(4,5))
for (idx in 1:ntraits){
  
  if (cnt > nplotspp){
    cnt = 1
    par(mfrow=c(4,5))
    print ("new page")
  }
  print(idx)
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  women <- merged_tab$gender_F1M2 == 1
  
  #if (colnames(merged_tab)[1] %in% cvd){
    if (length(unique(merged_tab[,1])) == 2){ # binary outcome
      print("binary")
      #merged_tab[,1] <- as.factor(merged_tab[,1])
      glm1 <- glm(merged_tab[,1] ~ merged_tab$age + merged_tab$gender_F1M2 + merged_tab$age*merged_tab$gender_F1M2, family=binomial(link="logit"))
      a <- anova(glm1, test = "Chisq")
      anova_p <- a[4,"Pr(>Chi)"]
      
      no_exremes <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
      glm2 <- glm(no_exremes[,1] ~ no_exremes$age + no_exremes$gender_F1M2 + no_exremes$age*no_exremes$gender_F1M2, family=binomial(link="logit"))
      a2 <- anova(glm2, test = "Chisq")
      anova_p2 <- a[4,"Pr(>Chi)"]
    } else{ # quantitative outcome
      a <-anova(lm(merged_tab[,1] ~ merged_tab$age + merged_tab$gender_F1M2 + merged_tab$age*merged_tab$gender_F1M2))
      anova_p <- a[3,"Pr(>F)"]
    
      no_exremes <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
      a2 <-anova(lm(no_exremes[,1] ~ no_exremes$age + no_exremes$gender_F1M2 + no_exremes$age*no_exremes$gender_F1M2))
      anova_p2 <- a2[3,"Pr(>F)"]
    }
  
  #plot only if p < 0.05
   #if (anova_p2 < 0.05 | anova_p < 0.05 ){
     #selected_plots <- append (selected_plots,plot_continuous(merged_tab, women, anova_p, anova_p2))
     plot_scatter_and_gam(merged_tab, women, anova_p, anova_p2)
     cnt <- cnt + 1
   #}

#}
}
dev.off()
ntraits



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

plot_scatter_and_gam <- function(merged_tab, women, anova_p, anova_p2){
  if (length(unique(merged_tab[,1])) == 2){ # binary outcome
    gammod <- gam(merged_tab[women,1]~s(merged_tab[women,"age"]), family=binomial(link = "logit"))
    gammod2 <- gam(merged_tab[!women,1]~s(merged_tab[!women,"age"]), family=binomial(link = "logit"))
  } else {
    gammod <- gam(merged_tab[women,1]~s(merged_tab[women,"age"]), method = "REML")
    gammod2 <- gam(merged_tab[!women,1]~s(merged_tab[!women,"age"]), method = "REML")
 }
  
  palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
  plot(merged_tab$age, merged_tab[,1],  col = merged_tab$gender_F1M2,  pch = 16, main = paste0(colnames(merged_tab)[1], "\np = ", format(anova_p, digits = 3),"\n20<age<75 p = ", format(anova_p2, digits = 3)), cex = 0.6, xlab = "age", ylab = colnames(merged_tab)[1])
  #points(merged_tab[merged_tab$smk == 1,"age"], merged_tab[merged_tab$smk == 1,1], pch = 1, col = col2transparent("gold", 125), cex = 0.8)
  palette(c("indianred1", "dodgerblue1"))
  
  lines(merged_tab[women,"age"],fitted(gammod),col="indianred1", lwd = 2)
  lines(merged_tab[!women,"age"],fitted(gammod2),col="royalblue1", lwd = 2)
  #p <- ggplot() +
  #  geom_point(aes(x=merged_tab$age, y=merged_tab[,1],color = as.factor(merged_tab$gender_F1M2), alpha = 0.1)) +
  #  geom_line(aes(x=merged_tab[women,"age"],y=fitted(gammod)), col="indianred1", size = 1.2) +
  #  geom_line(aes(x=merged_tab[!women,"age"],y=fitted(gammod2)), col="dodgerblue1", size = 1.2) +
  #  scale_color_manual(values=c("indianred1", "dodgerblue1")) +
  #  theme(legend.position="none", panel.background = element_rect(fill='white', colour='gray'), panel.grid.major = element_line(colour = "lightgray"))
  
  #return(p)
}




idx=11
merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
drugs <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_45drugs_1135metasubjects_atleast5users.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
drugs_m <- drugs[match(row.names(merged_tab), row.names(drugs), nomatch = 0),]
merged_tab_m <- merged_tab[match(row.names(drugs_m), row.names(merged_tab), nomatch = 0),]

women <- merged_tab_m$gender_F1M2 == 1

pdf("TFF3_drugs.pdf", width = 15, height = 15)
par(mfrow=c(4,5))
cnt = 1
for (d in 1:ncol(drugs_m))
{
  print(d)
  if (cnt > 20){
    cnt = 1
    par(mfrow=c(4,5))
    print ("new page")
  }
palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
plot(merged_tab_m$age, merged_tab_m[,1], col = merged_tab_m$gender_F1M2,  pch = 16, main = paste0(colnames(merged_tab_m)[1], " - ", colnames(drugs_m[d])), cex = 0.6, xlab = "age", ylab = colnames(merged_tab_m)[1])
palette(c("indianred1", "dodgerblue1"))
points(merged_tab_m[drugs_m[,d] == 1,"age"], merged_tab_m[drugs_m[,d] == 1,1], pch = 8, col = col2transparent("gold", 125), cex = 0.8)

cnt <- cnt + 1
}

dev.off()



palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
plot(merged_tab_m$age, merged_tab_m[,1], col = merged_tab_m$gender_F1M2,  pch = 16, main = paste0(colnames(merged_tab_m)[1], "\noral_contraceptive | anti_androgen_oral_contraceptive"), cex = 0.6, xlab = "age", ylab = colnames(merged_tab_m)[1])
palette(c("indianred1", "dodgerblue1"))
points( merged_tab_m[drugs_m[,"oral_contraceptive"] == 1 | drugs_m[,"anti_androgen_oral_contraceptive"] == 1,"age"], merged_tab_m[drugs_m[,"oral_contraceptive"] == 1 | drugs_m[,"anti_androgen_oral_contraceptive"] == 1,1], pch = 8, col = col2transparent("gold", 125), cex = 0.8)



