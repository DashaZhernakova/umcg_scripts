library(plyr)
library(ggplot2)
library(gridExtra)

args <- commandArgs(trailingOnly = TRUE)

source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
wd_path <- "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/"
setwd(wd_path)

# Phenotypes
traits_path <- args[1]
pheno_path <- args[3]
out_basepath <- args[3]
#tologtr <- c("LEU", "LY", "LYP", "MO", "MOP", "GR", "GRP", "BA", "BAP", "EO", "EOP", "ER", "TR", "TGL", "HAL1", "HALB", "AST", "ALT", "AF", "GGT", "LCRP", "TSH", "UKRO", "UKR24")
tologtr <- c()

traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
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

indices = 1:num_traits
cnt = 1
plots <- list()
res <- data.frame(matrix(nrow = num_traits, ncol = 10))
colnames(res) <- c("phenotype", "N", "N women", "N men", "mean (sd)", "mean (sd) women", "mean (sd) men", "median", "median women", "median men")

for (idx in indices){
  trait_id <- colnames(traits_m)[idx]
  print(trait_id)
  if (trait_id %in% tologtr){
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = T)
    trait_id <- paste0("log(", trait_id, "+1)")
  } else {
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = F)
  }
  merged_tab <- na.omit(merged_tab)
  merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
  colnames(merged_tab)[1] <- "phenotype"
  mu <- ddply(merged_tab, "gender_F1M2", summarise, grp.mean=mean(phenotype))
  stdev <- ddply(merged_tab, "gender_F1M2", summarise, grp.sd=sd(phenotype))
  med <- ddply(merged_tab, "gender_F1M2", summarise, grp.median=median(phenotype))
  p <- ggplot(merged_tab, aes(x = phenotype, fill = gender_F1M2, color = gender_F1M2)) + geom_density(alpha=0.4) + 
    scale_fill_manual(values=c("indianred1", "dodgerblue1")) +
    scale_color_manual(values=c("indianred1", "dodgerblue1")) +
    geom_vline(data=mu, aes(xintercept=grp.mean, color=gender_F1M2), linetype="dashed", size = 1) +
    theme_bw() + theme(legend.position="none") +
    xlab(trait_id)
  plots[[cnt]] <- p
  
  m <- formatC(mean(merged_tab$phenotype), digits = 3, format = "f")
  mw <- formatC(mu[1,"grp.mean"], digits = 3, format = "f")
  mm <- formatC(mu[2,"grp.mean"], digits = 3, format = "f")
  sdev <- formatC(sd(merged_tab$phenotype), digits = 3, format = "f")
  sdw <- formatC(stdev[1,"grp.sd"], digits = 3, format = "f")
  sdm <- formatC(stdev[2,"grp.sd"], digits = 3, format = "f")
  med_all <- formatC(median(merged_tab$phenotype), digits = 3, format = "f")
  medw <- formatC(med[1,"grp.median"], digits = 3, format = "f")
  medm <- formatC(med[2,"grp.median"], digits = 3, format = "f")
  
  res[cnt,] <- c(trait_id, nrow(merged_tab), nrow(merged_tab[merged_tab$gender_F1M2 == 1,]), nrow(merged_tab[merged_tab$gender_F1M2 == 2,]),
                 paste0(m, " (", sdev, ")") , paste0(mw, " (", sdw, ")"), paste0(mm, " (", sdm, ")"), med_all, medw, medm)
  cnt <- cnt + 1 
}
write.table(res, file = paste0(out_basepath, ".summary.txt"), sep = "\t", quote = F, col.names = NA)
ml <- marrangeGrob(grobs= lapply(plots, "+", theme(plot.margin=margin(20,20,20,20))), 
                   layout_matrix = matrix(1:20, 5, 4, TRUE), top=NULL)
ggsave(paste0(out_basepath, ".density.pdf"), ml, width = 20, height = 20)
  

