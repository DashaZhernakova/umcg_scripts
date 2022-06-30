setwd("/Users/dashazhernakova/Docuwoments/UMCG/data/gender_differences/omics/test_loess")
library("fANCOVA")
library(RColorBrewer)

#proeomics
#prot <- as.data.frame(t(read.table("/Users/dashazhernakova/Docuwoments/UMCG/data/olink/olink_shared_folder/Data/rawProteinData/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

#phenotypes
#prot <- read.table("../phenotypes_nonbinary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#diseases
#prot <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150715_diseases_1135patients.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#drugs
prot <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_45drugs_1135metasubjects_atleast5users.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#diet 
prot <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_diet_1135metasubjects_log_transform_imputed.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#pheno from Science
prot <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_intrinsicfactors_1135metasubjects.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#smoking
prot <- read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/phenotypeDataFromScience/20150618_smoke_1135metasubjects_noMissing.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

pheno <- as.data.frame(t(read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/correctedForAgeGenderSmokingContracCellCounts/age_gender_smk_contrac_cell_counts.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
match_idx = intersect(row.names(prot),row.names(pheno))
prot_m <- prot[match_idx,]
pheno_m <- pheno[match_idx,]
prot_m <- prot_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]
women <- pheno_m$gender_F1M2 == 1
nrow(pheno_m)
palette(brewer.pal(n = 8, name = "Set1"))

span_val = 0.75

nprot <- ncol(prot_m)
npages = 1
nplotspp <- ceiling(nprot/npages)
pdf("smoking_scatter_loess0.75.pdf", width = 15, height = 15)
for (page in 1:npages){
 par(mfrow=c(ceiling(nplotspp/5),5))

 for (np in 1:nplotspp){
   idx <- nplotspp*(page - 1) + np
   if (idx > nprot){
     break
   }
   print(idx)
   prot_na.rm <- prot_m[!is.na(prot_m[,idx]),idx]
   pheno_na.rm <- pheno_m[!is.na(prot_m[,idx]),]
   women <- pheno_na.rm$gender_F1M2 == 1
   
   l1 <- lm(prot_na.rm~pheno_na.rm$age)
   l2 <- lm(prot_na.rm~pheno_na.rm$age + pheno_na.rm$gender_F1M2)
   l3 <- lm(prot_na.rm~pheno_na.rm$age + pheno_na.rm$gender_F1M2 + pheno_na.rm$age*pheno_na.rm$gender_F1M2)
   a <- anova(l1,l2, l3)
   anova_p <- a[3,"Pr(>F)"]
   #main = paste0(colnames(prot_m)[idx], "; p=", formatC(anova_p, format = "e", digits = 2))
   plot(pheno_na.rm$age, prot_na.rm, col = pheno_na.rm$gender_F1M2,  pch = 16, main = colnames(prot_m)[idx], cex = 0.6, xlab = "age", ylab = colnames(prot_m)[idx])
   
   l <- loess(prot_na.rm[women]~pheno_na.rm[women,"age"], span = span_val)
   #l <- loess.as(pheno_na.rm[women,"age"], prot_na.rm[women], degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F)
   pl <- predict(l)
   
   lines(pl, x=pheno_na.rm[women,"age"], col = 1, lwd = 5)
   #gammod <- gam(prot_na.rm[women]~pheno_na.rm[women,"age"], family=binomial(link = "logit"))
   #lines(pheno_na.rm[women,"age"],fitted(gammod),col=1, lwd = 3)
   
   l <- loess(prot_na.rm[!women]~pheno_na.rm[!women,"age"], span = span_val)
   #l<-loess.as(pheno_na.rm[!women,"age"], prot_na.rm[women], degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F)
   pl <- predict(l)
   
   lines(pl, x=pheno_na.rm[!women,"age"], col = 2, lwd = 5)
   #gammod2 <- gam(prot_na.rm[!women]~pheno_na.rm[!women,"age"], family=binomial(link = "logit"))
   #lines(pheno_na.rm[!women,"age"],fitted(gammod2),col=2, lwd = 3)
   
 }
}
dev.off()


#prot <- as.data.frame(t(read.table("/Users/dashazhernakova/Docuwoments/UMCG/data/olink/LLDeep_expression/gene_level/LLD_genelevel_readcounts.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F))) 

#for (span_val in seq(0.1,1, length.out = 9)){
#  plot(pheno_na.rm$age, prot_na.rm, col = pheno_na.rm$gender_F1M2,  pch = 16, main = span_val, cex = 0.6, xlab = "age", ylab = colnames(prot_m)[idx])
#  
#  l <- loess(prot_na.rm[women]~pheno_na.rm[women,"age"], span = span_val)
#  #l <- loess.as(pheno_na.rm[women,"age"], prot_na.rm[women], degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F)
#  pl <- predict(l)
#  lines(pl, x=pheno_na.rm[women,"age"], col = 1, lwd = 5)
  
#  l <- loess(prot_na.rm[!women]~pheno_na.rm[!women,"age"], span = span_val)
#  #l<-loess.as(pheno_na.rm[!women,"age"], prot_na.rm[women], degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = F)
#  pl <- predict(l)
#  lines(pl, x=pheno_na.rm[!women,"age"], col = 2, lwd = 5)
  
#}
  