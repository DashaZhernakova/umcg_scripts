args <- commandArgs(trailingOnly = TRUE)
library(RColorBrewer)
library('dplyr')
library('mgcv')

source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")

######################################################################################
boxy = F

if (boxy){
  wd_path <- '/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/'
} else {
  wd_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/data/"
}

if  (length(args) > 1) {                                                                                                                                                                                                                       traits0_path <- args[1]
out_prefix <- args[2]
} else {
  #traits0_path <- "../LLD_expression/gene_read_counts_BIOS_and_LLD_passQC.only_LLD.no_zeros.TMM.CPM.tsv.gz.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"
  #traits0_path <- "LLD_expression_ageDEgenes.txt.gz"
  traits0_path <- "expression_selected.signif_inter.txt"
  
  out_prefix <- "test_ageDEgenes"
}
gte_path <- "gte_all.txt"
pheno_path <- "age_gender_smk_contrac_cell_counts.cc_outliers_na.txt"
gene_table_path <- "geneid_to_gene_proteincoding_mainchr.txt"
plot_basepath <- out_prefix
correct_for_cellcounts = TRUE
make_plots = TRUE

print(paste('wd:', wd_path))
print(paste('traits:',traits0_path))
print(paste('pheno:', pheno_path))
print(paste('out_prefix:', out_prefix))
print(paste('correct for cell counts:', correct_for_cellcounts))
print(paste('make plots:', make_plots))



setwd(wd_path)

# expression
traits0 <- as.data.frame(t(read.table(gzfile(traits0_path), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
gte <- read.table(gte_path, sep = "\t", as.is = T, check.names = F)
gte_m <- gte[match(row.names(traits0), gte[,2], nomatch = 0 ),]
traits <- traits0[match(gte_m[,2], row.names(traits0), nomatch = 0 ),]
all(row.names(traits) == gte_m[,2])
row.names(traits) = gte_m[,1]
#cvd <- c("relativesCVD","hr","p","p_axis","pq","qrs","qrs_axis","qt","qtc","t_axis","dbp","hbf","map","sbp","cho","crph","glu","hb1c","hbac","hdc","k","ldc","tgl","ins","homair","bmi","hip","angioplasty","arythmia","chestpain","diabetes","diabtype","dilataorta","edema","heartattack","heartfailure","highchol","hypertension","lof","narrowcarotis","stroke","T1Dy_n","T2Dy_n","RelativesDiab","added_compl_legs_1_yes_0_no","added_pain_hands_or_feet_1y_0n","added_stiff_hands_or_feet_1y_0n","feetwounds","jointpainfeet","jointpainhands","legcomplnight","legcomplwalking","stiffnessfeet","stiffnesshands")
#traits <- as.data.frame(t(read.table(gzfile("/Users/dashazhernakova/Documents/UMCG/data/olink/LLDeep_expression/gene_level/LLD_genelevel_readcounts.nozeros.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

num_traits <- ncol(traits)

pheno0 <- as.data.frame(t(read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
pheno <- na.omit(pheno0)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]



gene_table <- read.table(gene_table_path, header = T, sep = "\t", as.is = T, check.names = F)


#pheno_m$gender_F1M2 <- as.factor(pheno_m$gender_F1M2)

ntraits <- ncol(traits_m)
ntraits
nplotspp <-20


if (length(args) > 3) {
  start_idx <- args[3]
  end_idx <- args[4]
  plot_path <- paste0(plot_basepath, ".", start_idx, "-", end_idx,".pdf")
  
} else {
  start_idx <- 1
  end_idx <- ncol(traits_m)
  plot_path <- paste0(plot_basepath, ".pdf")
  
}

n_points = 300
res_dif_all <- data.frame(age = seq(20, 75, length = n_points))
res_summary <- data.frame()

cnt = 1
if (make_plots){
  pdf(plot_path, width = 15, height = 15)
  par(mfrow=c(5,4)) 
}
for (idx in start_idx:end_idx){
  if (cnt > nplotspp & make_plots){
    cnt = 1
    par(mfrow=c(5,4))
  }
  print(idx)
  
  trait_id <- colnames(traits_m)[idx]
  trait_name = gene_table[gene_table[,1] == trait_id,2]
  if (length(trait_name) > 0){ #if gene id in gene convertion table
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
    
    
    res_dif = NULL
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, correct_for_cellcounts, n_points, make_plots)
    
    res_summary[trait_id,'inter_p'] = res_dif_lst[["inter_p"]]
    res_summary[trait_id,'g_beta'] = res_dif_lst[["g_beta"]]
    res_summary[trait_id,'g_pv'] = res_dif_lst[["g_pv"]]
    
    if (res_dif_lst[["inter_p"]] < 0.05){
      #break
      cnt <- cnt + 1
      res_dif_all[,trait_id] <- res_dif_lst[["dif"]]
    }
  }
  
}
if (make_plots){
  dev.off()
}

write.table(res_dif_all, file = paste0(out_prefix, ".diff.txt"), sep = "\t", quote = F, col.names = NA)
write.table(res_summary, file = paste0(out_prefix, ".summary.txt"), sep = "\t", quote = F, col.names = NA)

