args <- commandArgs(trailingOnly = TRUE)
source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/factors")

rm_outliers <- function(merged_tab){
  merged_tab <- na.omit(merged_tab)
  w <- merged_tab[merged_tab$gender_F1M2 == 1,]
  m <- merged_tab[merged_tab$gender_F1M2 == 2,]
  mq1 <- quantile(m[,1], probs = 0.25)
  mq3 <- quantile(m[,1], probs = 0.75)
  miqr <- mq3 - mq1
  m_clean <- m[m[,1] < mq3 + 1.5*miqr & m[,1] > mq1 - 1.5*miqr,]
  
  wq1 <- quantile(w[,1], probs = 0.25)
  wq3 <- quantile(w[,1], probs = 0.75)
  wiqr <- wq3 - wq1
  w_clean <- w[w[,1] < wq3 + 1.5*wiqr & w[,1] > wq1 - 1.5*wiqr,]
  
  tab_nooutliers <- rbind(w_clean, m_clean)
  return(tab_nooutliers)
}


traits_path <- "../v4/data/LL_phenotypes_merged_all.log_some.v5.txt"
pheno_path <- "factors+med.txt"
pheno_name=args[1]
med=args[2]

cat("Data paths:\nphenotype traits:", traits_path, "\r\ncovariates:", pheno_path, "phenotype:", pheno_name, "medication name:", med, "\n")

# read phenotype traits of interest
traits <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

# read age, gender and other covariate phenotypes
pheno <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)


#order samples in the two tables
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)

if (! is.na(med)){
  cat("Medication will be used as a covariate\n")
  merged_tab <- cbind(traits_m[, pheno_name], pheno_m[,1:9], pheno_m[,med])
  colnames(merged_tab)[1] <- "phenotype"
  colnames(merged_tab)[ncol(merged_tab)] <- "med"
  
  merged_tab <- rm_outliers(merged_tab)
  
  preds <-  colnames(merged_tab)[c(4,5,8:ncol(merged_tab))]
  terms_binary <- " + SMK1 + SMK3 + s(age, by = SMK1) + s(age, by = SMK3)  + interaction(SMK1,gender_F1M2) + interaction(SMK3, gender_F1M2) + s(age, by = interaction(SMK1, gender_F1M2)) + s(age, by = interaction(SMK3, gender_F1M2)) + med + s(age, by = med) + interaction(med, gender_F1M2) + s(age, by = interaction(med, gender_F1M2))"
} else {
  cat("No medication to use\n")
  merged_tab <- cbind(traits_m[, pheno_name], pheno_m[,1:10])
  colnames(merged_tab)[1] <- "phenotype"

  merged_tab <- rm_outliers(merged_tab)
  
  preds <-  colnames(merged_tab)[c(4,5,8:ncol(merged_tab))]
  terms_binary <- " + SMK1 + SMK3 + s(age, by = SMK1) + s(age, by = SMK3)  + interaction(SMK1,gender_F1M2) + interaction(SMK3, gender_F1M2) + s(age, by = interaction(SMK1, gender_F1M2)) + s(age, by = interaction(SMK3, gender_F1M2))"

}
terms <-  paste0(" + s(", paste(preds, collapse = ")+ s("), ")")
terms_inter3 <- paste0(" + ti(", paste(preds, collapse = ", age, by = gender_F1M2)+ ti("), ", age, by = gender_F1M2)")
terms_inter_age <- paste0(" + s(", paste(preds, collapse = ", age)+ s("), ", age)")
terms_inter_sex <- paste0(" + s(", paste(preds, collapse = ", by = gender_F1M2)+ s("), ", by = gender_F1M2)")


full_formula <- as.formula(paste("phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2) ", terms, terms_inter3, terms_inter_age, terms_inter_sex, terms_binary, sep = " "))
print(full_formula)

full_fit <- gam(full_formula, data = merged_tab, method = "REML", select=T)

#fit <- gam(phenotype ~ SMK1 + SMK3 + gender_F1M2 + s(age) + s(LTE_SUM) + s(LDI_SUM) + s(total_mwk_VAL) + s(  total_scor_VAL) + s(MVPA_mwk_VAL) + s(MVPA_scor_VAL) +  s(LLDS_T1A) + s(SumOfalcohol),data = merged_tab, method = "REML", select=T)
s <- summary(full_fit)
write.table(s$p.table, file = paste0(out_path, pheno_name, ".p.table.txt"), sep = "\t", quote = F, col.names = NA)
write.table(s$s.table, file = paste0(out_path, pheno_name, ".s.table.txt"), sep = "\t", quote = F, col.names = NA)
#}


