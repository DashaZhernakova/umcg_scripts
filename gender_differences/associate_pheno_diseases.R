library('mgcv')


wd_path <- "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/"

# Phenotypes
traits_path <- "Laboratory_assessment_Blood_1A.dat"
st_col = 2

plot_basepath <- paste0("plots/", traits_path, ".pdf")

pheno_path <- "age_gender_all_LL.txt"
setwd(wd_path)

# traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits0 <- traits0[,seq(st_col,ncol(traits0))]
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)

# Age, gender and other phenotypes
pheno0 <- read.delim(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- na.omit(pheno0)

traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)
num_traits

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

# diseases
diseases <- read.delim("Questionnaire_Health_General_1A.selected.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits_m2 <- traits_m[match(row.names(diseases), row.names(traits_m), nomatch = 0),]
pheno_m2 <- pheno_m[match(row.names(diseases), row.names(pheno_m), nomatch = 0),]
diseases_m2 <- diseases[match(row.names(traits_m2), row.names(diseases), nomatch = 0),]

all(row.names(diseases_m2)==row.names(traits_m2))
all(row.names(diseases_m2)==row.names(pheno_m2))
all(row.names(traits_m2)==row.names(pheno_m2))


res_table_glm <- matrix(NA, ncol(diseases_m2), ncol(traits_m2))
res_table_gam <- matrix(NA, ncol(diseases_m2), ncol(traits_m2))

row.names(res_table_glm) <- row.names(res_table_gam) <- colnames(diseases_m2)
colnames(res_table_glm) <- colnames(res_table_gam) <- colnames(traits_m2)

for (i in 1:ncol(diseases_m2)){
  for (j in 1:ncol(traits_m2)){
    print (paste0(i, ",", j))
    merged_tab <- cbind(traits_m2[,j], diseases_m2[,i], pheno_m2)
    tr_name <- colnames(traits_m)[j]
    d_name <- colnames(diseases_m2)[i]
    colnames(merged_tab) <- c("phenotype", "disease", colnames(pheno_m2))
    row.names(merged_tab) <- row.names(pheno_m2)
    merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
    tryCatch({
    glm.fit <- glm(disease ~ phenotype + age + gender_F1M2, data = merged_tab, family = binomial(link = "logit"))
    glm_p <- summary(glm.fit)$coefficients["phenotype",4]
    res_table_glm[i,j] <- glm_p
    },error=function(e) {
          message(paste("Fitting failed for ", i, j))
     })
    #gam.fit <- gam(disease ~ gender_F1M2 + s(phenotype) + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML", family = binomial(link = "logit"))
    #gam_p <- summary(gam.fit)$s.pv[1]
    #res_table_gam[i,j] <- gam_p
  }
}

write.table(res_table_glm, file = "summary_tables/blood_pheno_vs_diseases_glm_quick.txt", sep = "\t", quote = F, col.names = NA)
#write.table(res_table_gam, file = "summary_tables/blood_pheno_vs_diseases_gam.txt"), sep = "\t", quote = F, col.names = NA)





