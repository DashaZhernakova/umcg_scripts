library(mgcv)
library(dplyr)

setwd("/groups/umcg-lifelines/tmp01/projects/ov20_0051/umcg-dzhernakova/gender_difs/prediction_model")

traits_path <- "../v5/data/LL_phenotypes_merged_all.log_some.v5.txt"
pheno_path <- "age_gender_new_CVD.v5.txt"

traits <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)
d <- cbind(pheno_m, traits_m)

d <- mutate(d, gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
d <- mutate(d, SMK3 = ordered(SMK3, levels = c('0', '1')))
d <- mutate(d, T2D = ordered(T2D, levels = c('0', '1')))

w <- d[d$gender_F1M2 == "1",]
m <- d[d$gender_F1M2 == "2",]


phenos <- colnames(traits_m)
sign_phenos <- c()
res <- data.frame(matrix(nrow = length(phenos), ncol = 7))
row.names(res) <- phenos
colnames(res) <- c("sex", "age", "pheno", "age*sex", "pheno*sex", "pheno*age", "pheno*age*sex")
for (p in phenos){
    d_cp <- d
    idx=which(colnames(d_cp) == p)
    colnames(d_cp)[idx] <- "pheno"
    if (length(unique(d_cp$pheno)) > 3){
        g <- gam(new_cvd~gender_F1M2 + s(age) + s(pheno) + s(age, by = gender_F1M2) + s(pheno, by = gender_F1M2) + ti(pheno, age) + ti(pheno, age, by = gender_F1M2), data = d_cp, method = "REML", family = binomial(link = "logit"), select = T)
        res[p,] <- c(summary(g)$p.table["gender_F1M2.L",4], summary(g)$s.table[,"p-value"])
    } else {
        g <- gam(new_cvd~gender_F1M2 + s(age) + pheno + s(age, by = gender_F1M2) + interaction(pheno, gender_F1M2) + s(age, by = pheno) + s(age, by = interaction(pheno, gender_F1M2)), data = d_cp, method = "REML", family = binomial(link = "logit"), select = T)
        res[p,] <- c(summary(g)$p.table["gender_F1M2.L",4], summary(g)$s.table["s(age)","p-value"], summary(g)$p.table["pheno.L",4], 
            summary(g)$s.table["s(age):gender_F1M22","p-value"], min(summary(g)$p.table[c(4,5,6),4], na.rm = T), summary(g)$s.table["s(age):pheno1","p-value"],
            min(summary(g)$s.table[c(4,5,6,7),"p-value"], na.rm = T))
    }
    print(summary(g))
 
}
print(res)

