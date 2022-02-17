library(rprojroot)
library(tidyverse)

  script_folder <- "C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/"
  config_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/config_local.yml"

cat("script folder:", script_folder, "\n")
source(paste0(script_folder, "/preprocessing_gam_fitting_functions.R"))
source(paste0(script_folder, "/get_breakpoints.R"))
source(paste0(script_folder, "/additional_functions.R"))
source(paste0(script_folder, "/plotting_functions.R"))

cat("Using config file: ", config_path, "\n")
config <- config::get(file = config_path)
# save the config in results folder

#
# Read data
#

traits_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/CVD3_olinkNormal_1447_LLDsamples_ProtNames.txt"
#traits_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/LLD_bloodlipids_nmr.txt"

pheno_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/LLD_covariates+phenotypes_from_full_LL.txt"
cat("Data paths:\nphenotype traits:", traits_path, "\r\ncovariates:", pheno_path, "\noutput base folder:", config$basedir_path, "\n\n")

# read phenotype traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)

traits <- na.omit(traits)
pca <- prcomp(traits, scale = T, center = T)


# read age, gender and other covariate phenotypes
pheno <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

#order samples in the two tables
pca_m <- pca$x[match(row.names(pheno), row.names(pca$x), nomatch = 0 ), ]
pheno_m <- pheno[match(row.names(pca_m), row.names(pheno), nomatch = 0), ]
all(row.names(pca_m) == row.names(pheno_m))
num_traits <- ncol(pca_m)
num_pheno <- ncol(pheno_m)

num_traits
num_pheno

res <- data.frame(matrix(ncol = 4, nrow = num_traits*num_pheno))
colnames(res) <- c("PC","pheno", "cor", "pval")
cnt <- 1
for (i in 1:num_traits){
  pc <- colnames(pca_m)[i]
  for (j in 1:num_pheno){
    p = colnames(pheno_m)[j]
    cur_pheno <- pheno_m[! is.na(pheno_m[,j]),]
    cur_pca <- pca_m[row.names(cur_pheno),i]
    cat(pc, p, "\n")
    res[cnt, 1] <- pc
    res[cnt, 2] <- p
    if (length(unique(cur_pheno[,j])) > 2){ # not binary
      c <- cor.test(cur_pca, cur_pheno[,j], method = "spearman")
      res[cnt, 3] <-c$estimate
      res[cnt, 4] <- c$p.value
    } else {
      c <- t.test(cur_pca ~ cur_pheno[,j])
      res[cnt, 3] <- NA
      res[cnt, 4] <- c$p.value
    }
    cnt <- cnt + 1
  }
}
res$padj <- p.adjust(res$pval, method = "BH")
write.table(res, "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/results/tables/Olink_PCA_vs_pheno.v2.txt", sep = "\t", quote = F, col.names = NA)



library(RColorBrewer)

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}
palette(c(col2transparent("#ff9999", 150),col2transparent("#99ccff", 150)))
par(tck = -.01, # Reduce tick length
    xaxs = "i", yaxs = "i")
plot(pca_m[,2], pca_m[,3], col = pheno_m[,"gender_F1M2"], xlab = "PC2", ylab = "PC3", 
     main = "sex is correlated with proteomics PC3", cex.main = 0.8,
     pch = 16, cex = 0.8, frame.plot = F, axes = T)
abline(h = pretty(pca_m[,3]), col = "grey90")
abline(v = pretty(pca_m[,2]), col = "grey90")
