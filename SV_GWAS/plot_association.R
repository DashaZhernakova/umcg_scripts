library(ggplot2)
library(patchwork)
library(plyr)
cohorts <- c("LLD", "500FG", "DAG3")
snp <- "9:136149830"
rsid <- "rs532436"
a1 <- "A"
a2 <- "G"

#snp="9:136155000"
#rsid="rs635634"
#a1 <- "T"
#a2 <- "C"
sv <- "F.prausnitzii:102"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:577_579"
svtype <- "dSV"
sv <- "F.prausnitzii:33"
sv_name <- "Faecalibacterium cf. prausnitzii KLE1255:885_887"
svtype <- "vSV"
c="DAG3"

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
    abo$ABO_O <- as.factor(abo$ABO_O )
    abo$Bloodtype <- as.factor(abo$Bloodtype )
    return(abo)
}

#pdf(paste0(snp, "-", sv, ".pdf"))
#par(mfrow=c(length(cohorts), 2))

#for (c in cohorts){
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
    
    abo <- read.delim(paste0("/data/umcg-tifn/SV/SV_GWAS/data/pheno/", c, "_Bloodtypes"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    abo <- format_abo_bloodgroup(abo)
    
    
    ids <- intersect(row.names(pheno), row.names(covar))
    
    if (c != "DAG3"){
        #m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "age", "read_number")], sex[ids,5], geno[ids,], fut2_geno[ids, "secr"], abo[ids, c("ABO_O", "Bloodtype")])
        #colnames(m) <- c("sv", "abundance", "PC1", "PC2", "age", "read_number", "sex", "geno_factor", "genotype", "FUT2_status", "ABO_O", "ABO_blood_group")
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "age", "read_number")], sex[ids,5], geno[ids,], fut2_geno[ids, "secr"])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2", "age", "read_number", "sex", "geno_factor", "genotype", "FUT2_status")


        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~  PC1 + PC2 + PC3 + age + read_number + sex")
        # SV vs ABO O group
        full_formula_ABO <- as.formula("sv ~ ABO_O + abundance + PC1 + PC2 + age + read_number + sex")

        # SV vs FUT2 status
        full_formula_FUT2 <- as.formula("sv ~ FUT2_status + abundance + PC1 + PC2 + age + read_number + sex")

        # SV vs ABO vs FUT2 
        full_formula_ABOFUT2 <- as.formula("sv ~ ABO_FUT2 + abundance + PC1 + PC2 + age + read_number + sex")
    } else {
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "PC3", "PC4", "PC5", "age", "read_number")], sex[ids,5], geno[ids,], fut2_geno[ids, "secr"], abo[ids, c("ABO_O", "Bloodtype")])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2",  "PC3", "PC4", "PC5", "age", "read_number", "sex", "geno_factor", "genotype", "FUT2_status", "ABO_O", "ABO_blood_group")
        
        
        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~ PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")

        # SV vs ABO O group
        full_formula_ABO <- as.formula("sv ~ ABO_O + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")

        # SV vs FUT2 status
        full_formula_FUT2 <- as.formula("sv ~ FUT2_status + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")

        # SV vs ABO vs FUT2 
        full_formula_ABOFUT2 <- as.formula("sv ~ ABO_FUT2 + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
    }
    m <- na.omit(m)
    m$ABO_FUT2 <- interaction(m$FUT2_status, m$ABO_blood_group)
    m$ABO_FUT2 <- factor(m$ABO_FUT2, levels = c("secretor.A", "secretor.B", "secretor.AB", "secretor.O", "non-secretor.A", "non-secretor.B", "non-secretor.AB", "non-secretor.O"))
    m$ABO_O_FUT2 <- interaction(m$FUT2_status, m$ABO_O)
    m$ABO_O_FUT2 <- factor(m$ABO_O_FUT2, levels = c("secretor.O", "secretor.A/B", "non-secretor.O", "non-secretor.A/B"))
    m$SNP_FUT2 <- interaction(m$FUT2_status, m$geno_factor)
    m$SNP_FUT2 <- factor(m$SNP_FUT2, levels = c(paste0("secretor", ".", levels(m$geno_factor)), paste0("non-secretor", ".", levels(m$geno_factor))))
    
    if (svtype == "dSV"){
        fit <- glm(full_formula, data = m, family = binomial(link = "logit"))
        fit_covar <- glm(covar_formula, data = m, family = binomial(link = "logit"))
        fit_covar_abund <- lm(covar_formula_abundance, data = m)

        print(table(m[,c("sv", "geno_factor")]))
        print(table(m[m$FUT2_status == "secretor", c("sv", "geno_factor")]))
        print(table(m[m$FUT2_status == "non-secretor", c("sv", "geno_factor")]))

        print(table(m[,c("sv", "ABO_blood_group")]))
        print(table(m[m$FUT2_status == "secretor", c("sv", "ABO_blood_group")]))
        print(table(m[m$FUT2_status == "non-secretor", c("sv", "ABO_blood_group")]))
        
    } else if (svtype == "vSV"){
        fit <- lm(full_formula, data = m, )
        fit_covar <- lm(covar_formula, data = m)
        fit_covar_abund <- lm(covar_formula_abundance, data = m)

        ddply(m, "SNP_FUT2", summarise, grp.mean=mean(sv))


    }
    print(summary(fit))
    

    #fit_ABO <- glm(full_formula_ABO, data = m, family = binomial(link = "logit"))
    #fit_FUT2 <- glm(full_formula_FUT2, data = m, family = binomial(link = "logit"))
    #fit_ABOFUT2 <- glm(full_formula_ABOFUT2, data = m, family = binomial(link = "logit"))

    m$sv_corrected <- residuals(fit_covar)
    
    
    
    
    m$abundance_corrected <- residuals(fit_covar_abund)

    p1 <- ggplot(m, aes(x = geno_factor, y = sv_corrected)) +
        geom_boxplot(color="#30455f", fill="#5f74e9", alpha=0.2) +
        labs(title=paste0(sv, " : ", rsid), x = "genotype", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none")

    p2 <- ggplot(m, aes(x = geno_factor, y = abundance_corrected)) +
        geom_boxplot(color="#30455f", fill="#5f74e9", alpha=0.2) +
        labs(title=paste0("F.prausnitzii corrected species abundance : ", rsid), x = "genotype", y = "F.prausnitzii abundance") +
        theme_classic() + theme(legend.position = "none")
    
    p3 <- ggplot(m, aes(x = ABO_blood_group, y = sv_corrected)) +
        geom_boxplot(color="#30455f", fill="#5f74e9", alpha=0.2) +
        labs(title=paste0(sv, " : ABO blood group"), x = "ABO blood group", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none")

    p4 <- ggplot(m, aes(x = FUT2_status, y = sv_corrected)) +
        geom_boxplot(color="#30455f", fill="#5f74e9", alpha=0.2) +
        labs(title=paste0(sv, " : FUT2 status"), x = "FUT2 status", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none")
    
    p5 <-ggplot(m, aes(x = ABO_FUT2, y = sv_corrected, color = "dodgerblue")) +
        geom_boxplot(color="#30455f", fill="#5f74e9", alpha=0.2) +
        labs(title=paste0(sv, " : ABO-FUT2"), x = "", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position = "none")

    p6 <-ggplot(m, aes(x = ABO_O_FUT2, y = sv_corrected, color = "dodgerblue")) +
        geom_boxplot(color="#30455f", fill="#5f74e9", alpha=0.2) +
        labs(title=paste0(sv, " : ABO O group vs FUT2"), x = "", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position = "none")
    
    p7 <-ggplot(m, aes(x = SNP_FUT2, y = sv_corrected, color = "dodgerblue")) +
        geom_boxplot(color="#30455f", fill="#5f74e9", alpha=0.2) +
        labs(title=paste0(sv, " : genotype vs FUT2"), x = "", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), legend.position = "none")
    
    pdf(paste0("ABO-", rsid, ".", sv, ".", svtype,".pdf"), width = 8, height = 19)
    par(mfrow=c(4, 2))

    (p1 + p2) / (p3 + p4) / (p5 + p6) / p7
    dev.off()
#}


if (sv_type == "dSV"){
    # add sdev
    m[m$sv == 1,"sv_text"] <- "deletion"
    m[m$sv == 2,"sv_text"] <- "no deletion"
    m$sv_text <- as.factor(m$sv_text)
    p1 <- ggplot(m, aes(x = geno_factor)) + 
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c("#d2d6db", "dodgerblue"), 0.5)) +
        theme_classic() + labs( x = rsid, y = paste0("fraction of samples with ", svtype))

    p2 <- ggplot(m, aes(x = ABO_blood_group)) + 
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c("#d2d6db", "dodgerblue"), 0.5)) +
        theme_classic() + labs(title =  "SV vs ABO blood group", y = paste0("fraction of samples with ", svtype))


    p3 <- ggplot(m, aes(x = SNP_FUT2)) + 
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c("#d2d6db", "dodgerblue"), 0.5)) +
        theme_classic() + labs(title = "SV vs genotype by FUT2 status", x = "", y = paste0("fraction of samples with ", svtype)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0))
    
    p4 <- ggplot(m, aes(x = ABO_O_FUT2)) + 
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c("#d2d6db", "dodgerblue"), 0.5)) +
        theme_classic() + labs(title =  "SV vs blood group by FUT2 status", x = "", y = paste0("fraction of samples with ", svtype)) +
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0))
    
    pdf(paste0("ABO-", rsid, ".", sv, ".", svtype,".pdf"), width = 8, height = 19)
    par(mfrow=c(4, 2))

    ps <- (p1 + p2) / p3 / p4
    ps + plot_annotation(title = paste0(svtype, " ", sv_name))
    dev.off()

}


 as.data.frame(table(m[,c("sv", "geno_factor", "FUT2_status")]))


 sec <- m[m$FUT2_status == "secretor",]
 non <- m[m$FUT2_status == "non-secretor",]
 summary(lm(sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = m))
 summary(lm(sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = sec))
 summary(lm(sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = non))

summary(lm(sv ~ ABO_blood_group + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = m))
 summary(lm(sv ~ ABO_blood_group + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = sec))
 summary(lm(sv ~ ABO_blood_group + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = non))

summary(lm(sv ~ ABO_O + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = m))
 summary(lm(sv ~ ABO_O + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = sec))
 summary(lm(sv ~ ABO_O + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = non))




sec <- m[m$FUT2_status == "secretor",]
 non <- m[m$FUT2_status == "non-secretor",]
 summary(glm(sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = m, family = binomial(link = "logit")))
 summary(glm(sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = sec, family = binomial(link = "logit")))
 summary(glm(sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = non, family = binomial(link = "logit")))

summary(glm(sv ~ ABO_blood_group + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = m, family = binomial(link = "logit")))
 summary(glm(sv ~ ABO_blood_group + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = sec, family = binomial(link = "logit")))
 summary(glm(sv ~ ABO_blood_group + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = non, family = binomial(link = "logit")))

summary(glm(sv ~ ABO_O + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = m, family = binomial(link = "logit")))
 summary(glm(sv ~ ABO_O + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = sec, family = binomial(link = "logit")))
 summary(glm(sv ~ ABO_O + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex, d = non, family = binomial(link = "logit")))













p1 <- ggplot(m, aes(x = geno_factor, y = sv, size = )) +
        geom_violin() +
        labs(title=paste0(sv, " : ", snp), x = "genotype", y = "corrected dSV frequency") +
        theme_classic() + theme(legend.position = "none")
        #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
    
    p2 <- ggplot(m, aes(x = geno_factor, y = sv, color = "dodgerblue")) +
        geom_boxplot() +
        labs(title=paste0(sv, " : ", snp), x = "genotype", y = "corrected dSV frequency") +
        theme_classic() + theme(legend.position = "none")
      

    p3 <- ggplot(m, aes(x = ABO_FUT2, y = sv, color = "dodgerblue")) +
        geom_violin() +
        labs(title=paste0(sv, " : ", snp), x = "ABO vs FUT2", y = "corrected dSV frequency") +
        theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")

    p4 <-ggplot(m, aes(x = ABO_FUT2, y = sv, color = "dodgerblue")) +
        geom_boxplot() +
        labs(title=paste0(sv, " : ", snp), x = "ABO vs FUT2", y = "corrected dSV frequency") +
        theme_classic() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "none")
    



c = "DAG3"
