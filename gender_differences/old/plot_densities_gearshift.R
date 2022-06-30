args <- commandArgs(trailingOnly = TRUE)

source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")

if (length(args) < 2){
  wd_path <- "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/"
} else {
  wd_path <- args[2]
}
# Phenotypes
#traits_path <- "Laboratory_assessment_Blood_1A.dat"
traits_path <- args[1]
st_col = 2

plot_basepath <- paste0("plots/densities_", traits_path, ".pdf")

pheno_path <- "age_gender_all_LL.txt"

correct_for_cellcounts = F
make_plots = T

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

# trait description
#pheno_code_path <- paste0("/groups/umcg-lifelines/tmp01/projects/phenotypes/metadata/", gsub("_1A.dat","_M.dat", traits_path))
#pheno_code_table <- read.table(pheno_code_path, header = T, sep = "\t", as.is = T, check.names = F)

nplotspp = 20
n_points = 300
res_dif_all <- data.frame(age = seq(20, 75, length = n_points))
res_summary <- data.frame()

cnt = 1
pdf(plot_basepath, width = 15, height = 15)
par(mfrow=c(5,4)) 

indices = 1:ncol(traits_m)
cnt = 1

for (idx in indices){
  if (cnt > nplotspp & make_plots){
    cnt = 1
    par(mfrow=c(5,4))
  }
  print(idx)
  
  trait_id <- colnames(traits_m)[idx]
  trait_name = trait_id
  trait_descr <- pheno_code_table[pheno_code_table[,1] == trait_id,3]

  if (length(trait_name) > 0 & length(unique(traits_m[,idx])) > 1){ #if gene id in gene convertion table
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = F)
    plot(density(merged_tab[,1]), main = trait_name)
  }
  
}

dev.off()
