pdf("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/v3/results/plots/pregnancies.pdf", width = 10, height = 10)
par(mfrow=c(3,3))
indices = 1:ncol(traits_m)
for (idx in indices){
  trait_name <- colnames(traits_m)[idx]
  cat(idx, " : ", trait_name, "\n")
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
  colnames(merged_tab)[1] <- "phenotype"
  postmeno <- merged_tab[merged_tab$age > 55,]
  lmfit <- lm(phenotype ~ FEM13, data = postmeno)
  pv <- summary(lmfit)$coefficients["FEM13",4]
  boxplot(phenotype ~ FEM13, data = postmeno, main = paste0(trait_name, "\n", format(pv, digits = 3)))
}

dev.off()