library(ggplot2)
library(patchwork)
library(plyr)


cohorts <- c("LLD", "500FG", "DAG3")
snp <- "9:136141870"
rsid <- "rs2519093"
a1 <- "T"
a2 <- "C"

sv <- "F.prausnitzii:102"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:577_579"
svtype <- "dSV"

fut2_snp <- "19:49206674"
recode_genotypes <- function(geno, a1, a2){
    geno$ac <- geno[,1]
    geno$ac <- sub(paste0(a2, "/", a2), 0, geno$ac)
    geno$ac <- sub(paste0(a2, "/", a1), 1, geno$ac)
    geno$ac <- sub(paste0(a1, "/", a2), 1, geno$ac)
    geno$ac <- sub(paste0(a1, "/", a1), 2, geno$ac)
    geno$ac <- as.numeric(geno$ac)
    geno[,1] <- as.factor(geno[,1])
    return(geno)
}

get_fut2_secretor_status <- function(fut2_geno){
    fut2_geno$secr <- "secretor"
    fut2_geno$secr[fut2_geno[,"19:49206674"] == "A/A"] <- "non-secretor"
    fut2_geno$secr <- as.factor(fut2_geno$secr)
    return (fut2_geno)
}

format_abo_bloodgroup <- function(abo){
    abo$ABO_O <- "A/B"
    abo[abo$Bloodtype == "O", "ABO_O"] <- "O"
    abo$ABO_O <- as.factor(abo$ABO_O)
    abo$Bloodtype <- as.factor(abo$Bloodtype)
    abo$ABO_A <- "O/B"
    abo[abo$Bloodtype %in% c("AB", "A"), "ABO_A"] <- "A"
    abo$A_geno <- 0
    abo[abo$Blood_genotype %in% c("A", "AB"), "A_geno"] <- 1
    abo[abo$Blood_genotype == "AA", "A_geno"] <- 2
    abo$ABO_blood_group <- abo$Bloodtype
    return(abo)
}



for (c in cohorts){
    print(c)
    pheno <- read.delim(paste0("/data/umcg-tifn/SV/SV_GWAS/data/", c, ".", svtype, ".filtered.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    covar <- read.delim(paste0("/data/umcg-tifn/SV/SV_GWAS/data/", c, ".covariates.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    if (svtype == "dSV") pheno[,sv] <- as.factor(pheno[,sv])
    geno <- as.data.frame(t(read.delim(paste0("/data/umcg-tifn/SV/SV_GWAS/genotypes/", c, "/text_genotypes/", c, ".", snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
    geno <- recode_genotypes(geno, a1, a2)
    fut2_geno <- as.data.frame(t(read.delim(paste0("/data/umcg-tifn/SV/SV_GWAS/genotypes/", c, "/text_genotypes/", c, ".", fut2_snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
    fut2_geno <- get_fut2_secretor_status(fut2_geno)

    sex <- read.table(paste0("/data/umcg-tifn/SV/SV_GWAS/genotypes/", c, "/", c , "_filtered.fam"), header = F,  as.is = T, check.names = F)
    row.names(sex) <- sex[,2]

    abo <- read.delim(paste0("/data/umcg-tifn/SV/SV_GWAS/data/pheno/", c, ".abo_blood_group.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    abo <- format_abo_bloodgroup(abo)

    ids <- intersect(row.names(pheno), row.names(covar))

    if (c != "DAG3"){
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "age", "read_number")], sex[ids,5], geno[ids,], fut2_geno[ids, "secr"], abo[ids, ])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2", "age", "read_number", "sex", "geno_factor", "genotype", "FUT2_status", colnames(abo))

        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~  PC1 + PC2  + age + read_number + sex")
    } else {
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "PC3", "PC4", "PC5", "age", "read_number")], sex[ids,5], geno[ids,], fut2_geno[ids, "secr"], abo[ids, ])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2",  "PC3", "PC4", "PC5", "age", "read_number", "sex", "geno_factor", "genotype", "FUT2_status", colnames(abo))

        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~ PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
    }
    m <- na.omit(m)
