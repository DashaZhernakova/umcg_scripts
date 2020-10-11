library('mgcv')
args <- commandArgs(trailingOnly = TRUE)

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
diseases_path=args[1]
diseases <- read.delim(diseases_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
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

res_table_long <- matrix(NA, nrow = ncol(diseases_m2)*ncol(traits_m2), ncol = 8)
colnames(res_table_long) <- c("pheno1", "pheno2", "glm_beta", "glm_pval", "glm2_beta", "glm2_pval", "gam_beta", "gam_pval")
cnt = 1
for (i in 2:ncol(diseases_m2)){
  for (j in 1:ncol(traits_m2)){
    print (paste0(i, ",", j))
    merged_tab <- cbind(traits_m2[,j], diseases_m2[,i], pheno_m2)
    tr_name <- colnames(traits_m)[j]
    d_name <- colnames(diseases_m2)[i]
    colnames(merged_tab) <- c("phenotype", "disease", colnames(pheno_m2))
    row.names(merged_tab) <- row.names(pheno_m2)
    merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
    merged_tab$disease <- as.factor(merged_tab$disease)
    if (nrow(merged_tab[merged_tab$disease != 0,]) > 50){

    tryCatch({
    glm.fit <- glm(disease ~ phenotype + age + gender_F1M2, data = merged_tab, family = binomial(link = "logit"))
    glm_p <- summary(glm.fit)$coefficients["phenotype",4]
    glm_b <- summary(glm.fit)$coefficients["phenotype",1]
    
    glm.fit2 <- glm(disease ~ phenotype + age + gender_F1M2 + age*gender_F1M2, data = merged_tab, family = binomial(link = "logit"))
    glm_p2 <- summary(glm.fit2)$coefficients["phenotype",4]
    glm_b2 <- summary(glm.fit2)$coefficients["phenotype",1]
    
    if (summary(glm.fit2)$coefficients["age:gender_F1M22", 4] < 0.05) print(paste(tr_name, d_name))
    
    gam.fit <- gam(disease ~ gender_F1M2 + phenotype + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML", family = binomial(link = "logit"))
    gam_p <- summary(gam.fit)$p.pv["phenotype"]
    gam_b <- gam.fit$coefficients["phenotype"]
    
    },error=function(e) {
          message(paste("Fitting failed for ", i, j))
     })
    #gam.fit <- gam(disease ~ gender_F1M2 + s(phenotype) + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML", family = binomial(link = "logit"))
    #gam_p <- summary(gam.fit)$s.pv[1]
    #res_table_gam[i,j] <- gam_p
    
    res_table_long[cnt, "pheno1"] <- d_name
    res_table_long[cnt, "pheno2"] <- tr_name
    res_table_long[cnt, "glm_beta"] <- glm_b
    res_table_long[cnt, "glm_pval"] <- glm_p
    res_table_long[cnt, "glm2_beta"] <- glm_b2
    res_table_long[cnt, "glm2_pval"] <- glm_p2
    res_table_long[cnt, "gam_beta"] <- gam_b
    res_table_long[cnt, "gam_pval"] <- gam_p
    cnt = cnt + 1
    
   } 
  }
}
dis_name <- basename(diseases_path)
write.table(res_table_long, file = paste0("summary_tables/blood_pheno_vs_diseases_association2.",dis_name,".txt"), sep = "\t", quote = F, col.names = NA)
#write.table(res_table_gam, file = "summary_tables/blood_pheno_vs_diseases_gam.txt"), sep = "\t", quote = F, col.names = NA)






