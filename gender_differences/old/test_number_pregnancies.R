setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs")
traits_path <- "LL_phenotypes_merged_all.log_some.v2.txt"
pheno_path <- "pheno_tables/age_sex_pregnancies_menopause_G03fem.txt"
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)

pheno0 <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno0 <- pheno0[,c("age", "gender_F1M2", "FEM13")] 
pheno <- na.omit(pheno0)

#order samples in the two tables
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)

#traits_m <- traits_m[order(pheno_m$age), , drop = F]
#pheno_m <- pheno_m[order(pheno_m$age),]

cat("Number of available phenotypes: ", num_traits, "\n")
cat("Number of shared samples: ", nrow(traits_m), "\n")



pheno_m[pheno_m$FEM13 > 5, "FEM13"] <- 5

pdf("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v3/results/plots/pregnancies_corrected_for_age_sq_resid.pdf", width = 10, height = 10)
par(mfrow=c(3,3))
indices = 1:ncol(traits_m)
for (idx in indices){
  trait_name <- colnames(traits_m)[idx]
  cat(idx, " : ", trait_name, "\n")
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "zscore", log_tr = F, scale_tr = F)
  colnames(merged_tab)[1] <- "phenotype"
  postmeno <- merged_tab[merged_tab$age > 60,]
  lmfit <- lm(phenotype ~ FEM13 + age, data = postmeno)
  pheno_resid <- residuals(lmfit)
  pv <- summary(lmfit)$coefficients["FEM13",4]
  boxplot(pheno_resid ~ postmeno$FEM13, main = paste0(trait_name, "\n", format(pv, digits = 3)))
}

dev.off()
