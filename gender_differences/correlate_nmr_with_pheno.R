library(reshape2)
library(corrplot)
setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/")
traits_path <- "results/data/LLD_bloodlipids_nmr.txt"
pheno_path <- "results/data/LLD_covariates_080321.txt"
diseases_path <- "data/LLD_pheno_selection.290721.txt"
diseases_path <- "results/data/disease_selection.txt"
diseases_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/LLD_diseases_1660.txt"
diseases_path <- "C:/Users/Dasha/work/UMCG/data/MR/results2/AA_T2D/data/LLD/traits/int_raw_no_outlier3sd.txt"
dis_name <- ""
binary_pheno <- FALSE
binary_pheno <- TRUE
remove_outliers <- function(m){
  m <- na.omit(m)
  q1 <- quantile(m$phenotype, probs = 0.25)
  q3 <- quantile(m$phenotype, probs = 0.75)
  iqr <- q3 - q1
  m_clean1 <- m[m$phenotype < q3 + 1.5*iqr & m$phenotype > q1 - 1.5*iqr,]
  
  q12 <- quantile(m_clean1$disease, probs = 0.25)
  q32 <- quantile(m_clean1$disease, probs = 0.75)
  iqr2 <- q32 - q12
  m_clean2 <- m_clean1[m_clean1$disease < q32 + 1.5*iqr2 & m_clean1$disease > q12 - 1.5*iqr2,]
}


traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)

# Age, gender and other phenotypes
pheno0 <- read.delim(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- pheno0[,c("age", "gender_F1M2", "bmi")]
pheno <- na.omit(pheno)

traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)
num_traits

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

# diseases
diseases <- read.delim(diseases_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
if (! binary_pheno){
  rows <- row.names(diseases)
  diseases <- sapply(diseases, function(x) as.numeric(as.character(x)))
  row.names(diseases) <- rows
}

traits_m2 <- traits_m[match(row.names(diseases), row.names(traits_m), nomatch = 0),]
pheno_m2 <- pheno_m[match(row.names(diseases), row.names(pheno_m), nomatch = 0),]
diseases_m2 <- diseases[match(row.names(traits_m2), row.names(diseases), nomatch = 0),]

all(row.names(diseases_m2)==row.names(traits_m2))
all(row.names(diseases_m2)==row.names(pheno_m2))
all(row.names(traits_m2)==row.names(pheno_m2))

#for (d in c("IL.1.Beta","IL.6","IL.8","IL.10","IL.12P70","TNF.Alpha")){
#  diseases_m2[,d] <- log(diseases_m2[,d])
#}

#col2excl <- c("IL.6", "TNF.Alpha", "Biochem_HDL", "Biochem_LDL", "Biochem_TG", "Biochem_Cholesterol")
#diseases_m2 <- diseases_m2[ , -which(colnames(diseases_m2) %in% col2excl)]


####
##
####

res_table_long <- as.data.frame(matrix(NA, nrow = ncol(diseases_m2)*ncol(traits_m2), ncol = 5), stringsAsFactors = F)
colnames(res_table_long) <- c("pheno1", "pheno2", "beta", "pval", "qval")
cnt = 1
for (i in 1:ncol(diseases_m2)){
  #if(nrow(diseases_m2[!is.na(diseases_m2[,i]),]) > 500){
    for (j in 1:ncol(traits_m2)){
      print (paste0(i, ",", j))
      merged_tab <- cbind(traits_m2[,j], diseases_m2[,i], pheno_m2)
      tr_name <- colnames(traits_m)[j]
      d_name <- colnames(diseases_m2)[i]
      colnames(merged_tab) <- c("phenotype", "disease", colnames(pheno_m2))
      row.names(merged_tab) <- row.names(pheno_m2)
      merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
      merged_tab$phenotype <- scale(merged_tab$phenotype)
      if (binary_pheno){
        merged_tab$disease <- as.factor(merged_tab$disease)
        
        glm.fit <- glm(disease ~ phenotype + age + gender_F1M2 + bmi, data = merged_tab, family = binomial(link = "logit"))
        p <- summary(glm.fit)$coefficients["phenotype",4]
        b <- summary(glm.fit)$coefficients["phenotype",1]
      } else {
        merged_tab <- remove_outliers(merged_tab)
        
        merged_tab$disease <- scale(merged_tab$disease)
        lm.fit <- lm(disease ~ phenotype + age + gender_F1M2 + bmi, data = merged_tab)
        p <- summary(lm.fit)$coefficients["phenotype",4]
        b <- summary(lm.fit)$coefficients["phenotype",1]
      }  
          
      res_table_long[cnt, "pheno1"] <- d_name
      res_table_long[cnt, "pheno2"] <- tr_name
      res_table_long[cnt, "beta"] <- b
      res_table_long[cnt, "pval"] <- p
      cnt = cnt + 1
        
    } 
  #}
}

res_table_long[,"qval"] <- p.adjust(res_table_long[,"pval"], method = "BH")
res_table_long <- na.omit(res_table_long)

res_table_long <- res_table_long[res_table_long$pval < 0.05,]
write.table(res_table_long, file = paste0("results/results/tables//pheno_vs_nmr_association2.",dis_name,".txt"), sep = "\t", quote = F, col.names = NA)
#write.table(res_table_gam, file = "summary_tables/blood_pheno_vs_diseases_gam.txt"), sep = "\t", quote = F, col.names = NA)


res_table_wide <- dcast(res_table_long, pheno2 ~ pheno1, value.var = "beta")
res_table_pval <- dcast(res_table_long, pheno2 ~ pheno1, value.var = "qval")
row.names(res_table_wide) <- res_table_wide$pheno2
res_table_wide <- res_table_wide[,-1]
row.names(res_table_pval) <- res_table_pval$pheno2
res_table_pval <- res_table_pval[,-1]

res_table_wide <- res_table_wide[,colnames(diseases_m2)]
res_table_pval <- res_table_pval[, colnames(diseases_m2)]

clusters <- read.delim("results/results/tables/nmr_clustering_kmeans.k5.txt", as.is = T, check.names = F)
clusters <- clusters[order(clusters$cl),]
lipid_order <- clusters[clusters$signif_inter == T,1]
lipid_names <- paste0(clusters[clusters$signif_inter == T,"cl"], ":", clusters[clusters$signif_inter == T,1])

pdf("corrplot_nmr_pheno2.pdf", width = 5, height = 10)
m <- res_table_wide[lipid_order,]
row.names(m) <- lipid_names
corrplot(as.matrix(m), is.corr = F, 
         tl.col = "black", tl.cex = 0.8, insig = "label_sig", pch = 1, pch.cex = 0.5,
         p.mat = as.matrix(res_table_pval[lipid_order,]), na.label = " ",
         sig.level = 0.05)
dev.off()

pdf("corrplot_nmr_diseases2.pdf", width = 5, height = 10)
m <- res_table_wide[lipid_order,]
row.names(m) <- lipid_names
corrplot(as.matrix(m), is.corr = F,
         tl.col = "black", tl.cex = 0.7,insig = "label_sig", pch = 1, pch.cex = 0.5,
         p.mat = as.matrix(res_table_pval[lipid_order,]), 
         sig.level = 0.05)
dev.off()





setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/")
traits_path <- "C:/Users/Dasha/work/UMCG/data/MR/results2/AA_T2D/data/300OB/traits/nmr_raw_noOuterlier.txt"
pheno_path <- "C:/Users/Dasha/work/UMCG/data/MR/results2/AA_T2D/data/300OB/traits/demo_raw_noOuterlier.txt"
diseases_path <- "C:/Users/Dasha/work/UMCG/data/MR/results2/AA_T2D/data/300OB/traits/cvd_raw_noOuterlier.txt"

pheno <- pheno0[,c("age", "sex", "BMI")]
colnames(pheno) <- c("age", "gender_F1M2", "bmi")

cols_log = c(3,5,7,9,14,17,18,19,20)
cols_int <- c(10,11,12)

for (i in cols_log){
  diseases_m2[,i] <- log(diseases_m2[,i])
}
# (i in cols_int){
#  
#  diseases_m2[,i] <- qnorm((rank(diseases_m2[,i],na.last="keep")-0.5)/sum(!is.na(diseases_m2[,i])))
#}
diseases_m2 <- diseases_m2[ ,-cols_int]

pdf("corrplot_nmr_300ob.pdf", width = 5, height = 10)
m <- res_table_wide[lipid_order,]
row.names(m) <- lipid_names
corrplot(as.matrix(m), is.corr = F, 
         tl.col = "black", tl.cex = 0.8, insig = "label_sig", pch = 1, pch.cex = 0.5,
         p.mat = as.matrix(res_table_pval[lipid_order,]), na.label = " ",
         sig.level = 0.05)
dev.off()

