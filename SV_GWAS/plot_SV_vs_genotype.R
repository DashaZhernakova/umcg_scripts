library(ggplot2)
library(patchwork)
library(plyr)


cohorts <- c("LLD", "500FG", "DAG3")
snp <- "12:270805"
rsid <- "rs3817017"
a1 <- "T"
a2 <- "C"
sv <- "A.putredinis:10"
sv_name <- "Alistipes putredinis DSM 17216:133_134"
svtype <- "dSV"


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

make_plots<- function(m, svtype, sv_name, pval){
    if (svtype == "dSV"){
    colors = c("#d2d6db", "#225EA8", "#CC4C02")
    p1 <- ggplot(m, aes(x = geno_factor)) +
        geom_bar(aes(fill = sv_text), position = "fill") +
        scale_fill_manual(sv, values = alpha( c(colors), 0.5)) +
        theme_classic() + 
        labs(title=paste0(c, " ", sv, " : ", rsid, "\nP = ", formatC(pval, digits = 4)), x = rsid, y = paste0("fraction of samples with ", svtype))
	} else if (svtype == "vSV"){
	colors = c("#225EA8", "#CC4C02")

    p1 <- ggplot(m, aes(x = geno_factor, y = sv_corrected)) +
        geom_violin(trim=FALSE, color = colors[1]) +
        geom_boxplot(width=0.5, color = colors[1]) +
        labs(title=paste0(c, " ", sv, " : ", rsid, "\nP = ", formatC(pval, digits = 4)), x = "genotype", y = paste0(svtype, " corrected abundance")) +
        theme_classic() + theme(legend.position = "none") +
        scale_color_manual(values = c(colors[1]))
	}
    
    return(p1)

}


print_association_p <- function(m, svtype){

    covars = "abundance + PC1 + PC2 + age + read_number + sex"
    if (c == "DAG3"){
        covars = "abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex"
    }

    if (svtype == "vSV"){
        lm_g <- lm(as.formula(paste0("sv ~ genotype + ", covars)), d = m)
        
    } else if (svtype == "dSV"){
        lm_g <- glm(as.formula(paste0("sv ~ genotype + ", covars)), d = m, family = binomial(link = "logit"))
    }
    return(summary(lm_g)$coefficients["genotype",4])
}


plot_list <- list()
cnt <- 1
for (c in cohorts){
    print(c)
    pheno <- read.delim(paste0("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/data/", c, ".", svtype, ".filtered.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    covar <- read.delim(paste0("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/data/", c, ".covariates.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
    if (svtype == "dSV") pheno[,sv] <- as.factor(pheno[,sv])
    geno <- as.data.frame(t(read.delim(paste0("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/genotypes/", c, "/text_genotypes/", c, ".", snp, ".genotypes.txt"), header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)))
    geno <- recode_genotypes(geno, a1, a2)
    
    sex <- read.table(paste0("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/umcg-dzhernakova/SV_GWAS/genotypes/", c, "/", c , "_filtered.fam"), header = F,  as.is = T, check.names = F)
    row.names(sex) <- sex[,2]
    ids <- intersect(row.names(pheno), row.names(covar))

    if (c != "DAG3"){
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "age", "read_number")], sex[ids,5], geno[ids,])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2", "age", "read_number", "sex", "geno_factor", "genotype")

        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~  PC1 + PC2  + age + read_number + sex")
    } else {
        m <- cbind(pheno[ids,sv], covar[ids, c(gsub(":[0-9]+","", sv), "PC1", "PC2", "PC3", "PC4", "PC5", "age", "read_number")], sex[ids,5], geno[ids,])
        colnames(m) <- c("sv", "abundance", "PC1", "PC2",  "PC3", "PC4", "PC5", "age", "read_number", "sex", "geno_factor", "genotype")

        # SV vs genotype
        full_formula <- as.formula("sv ~ genotype + abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula <- as.formula("sv ~ abundance + PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
        covar_formula_abundance <- as.formula("abundance ~ PC1 + PC2 + PC3 + PC4 + PC5 + age + read_number + sex")
    }
    m <- na.omit(m)

    
    if (svtype == "dSV"){
        fit <- glm(full_formula, data = m, family = binomial(link = "logit"))
        fit_covar <- glm(covar_formula, data = m, family = binomial(link = "logit"))
        fit_covar_abund <- lm(covar_formula_abundance, data = m)
        

    } else if (svtype == "vSV"){
        fit <- lm(full_formula, data = m, )
        fit_covar <- lm(covar_formula, data = m)
        fit_covar_abund <- lm(covar_formula_abundance, data = m)
    }

    m$sv_corrected <- residuals(fit_covar)
    m$abundance_corrected <- residuals(fit_covar_abund)
	pval <- print_association_p(m, svtype)
	
	if (svtype == "dSV"){
		# add sdev
		m[m$sv == 1, "sv_text"] <- "deletion"
		m[m$sv == 2, "sv_text"] <- "no deletion"
		m$sv_text <- as.factor(m$sv_text)	
	}
	plot_list[[cnt]] <- make_plots(m, svtype, sv_name, pval)
	cnt <- cnt + 1
}


pdf(paste0("all_cohorts.", rsid, ".", sv, ".", svtype,".pdf"), width = 15, height = 5)
plot_list[[1]] + plot_list[[2]] + plot_list[[3]]
dev.off()