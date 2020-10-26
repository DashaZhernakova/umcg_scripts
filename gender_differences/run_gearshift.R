args <- commandArgs(trailingOnly = TRUE)

source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/get_breakpoints.R")
source("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/scripts/umcg_scripts/gender_differences/calc_derivatives.R")
wd_path <- "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/"

# Phenotypes
traits_path <- args[1]
out_basepath <- args[2]
st_col = 1

#pheno_path <- "age_gender_all_LL.txt"
pheno_path <- "age_gender_bmi_smk_all_LL.txt"

correct_for_cellcounts = F
make_plots = T
add_breakpoints = F

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

nplotspp = 20
n_points = 300
min_age = 20
max_age = 80
ttest_cutoff <- 4
deriv_cutoff <- 0.0004
covariateslinear <- c("SMK3")
covariatesnonlinear <- c("BMI")
#covariateslinear <- c()
#covariatesnonlinear <- c()
#out_basepath <- basename(traits_path)
#plot_basepath <- paste0("plots/", out_basepath, "breakpoints_intervals_t",ttest_cutoff,"_d", deriv_cutoff,".png")
plot_basepath <- paste0("plots/", basename(out_basepath), ".png")

res_dif_all <- data.frame(age = seq(min_age, max_age, length = n_points))
res_pred_all <- data.frame(age = c(seq(min_age, max_age, length = n_points), seq(min_age, max_age, length = n_points)))
res_summary <- data.frame()

cnt = 1
if (make_plots){
  nrows <- ceiling(sqrt(num_traits))
  size <- 3*nrows
  png(plot_basepath, width = size, height = size, units = 'in', res = 400)
  par(mfrow=c(nrows,nrows))
}


indices = 1:ncol(traits_m)
cnt = 1

for (idx in indices){
  trait_id <- colnames(traits_m)[idx]
  trait_name = trait_id
  print(idx)
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "zscore", log_tr = F)
  print(paste0(trait_name, ",", nrow(merged_tab)))
  res_dif = NULL
  res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name,  covariates_linear = covariateslinear, covariates_nonlinear = covariatesnonlinear, n_points = n_points, make_plots = make_plots, gam_family = gaussian(), label = '', min_age = min_age, max_age = max_age, add_breakpoints = add_breakpoints, t_threshold = ttest_cutoff, derivatives_cutoff = deriv_cutoff)
   
  if (res_dif_lst[["inter_p"]] < 0.05){
      #cnt <- cnt + 1
      res_dif_all[,trait_id] <- res_dif_lst[["dif"]]
      res_pred_all[,trait_id] <- res_dif_lst[["pdat"]]$pred
  }
  res_summary[trait_id,'inter_p'] = res_dif_lst[["inter_p"]]
  res_summary[trait_id,'g_beta'] = res_dif_lst[["g_beta"]]
  res_summary[trait_id,'g_pv'] = res_dif_lst[["g_pv"]]
  sex_dif_pval <- calculate_sex_diff_ttest(merged_tab, covariates = c(covariateslinear, covariatesnonlinear), min_age, max_age)
  res_summary[trait_id,'g_ttest_pv'] = sex_dif_pval
  res_summary[trait_id,'cohen_f2'] = res_dif_lst[["cohen_f2"]]
  
  if (!is.null(res_dif_lst[['breakpoints_intervals']])){
      res_summary[trait_id, "breakpoints_men"] = ifelse(length(res_dif_lst[['breakpoints_intervals']][[2]]) > 0, paste0(res_dif_lst[['breakpoints_intervals']][[2]], collapse = ","), "NA")
      res_summary[trait_id, "breakpoints_women"] = ifelse(length(res_dif_lst[['breakpoints_intervals']][[1]]) > 0, paste0(res_dif_lst[['breakpoints']][[1]], collapse = ","), "NA")
   }  
}
#write.table(res_dif_all, file=paste0("summary_tables/", traits_path, "diff.txt"), sep = "\t", quote = F, col.names = NA)
#write.table(res_pred_all, file=paste0("summary_tables/", traits_path, ".fitted.txt"), sep = "\t", quote = F, col.names = NA)
res_summary$inter_p_adj <- p.adjust(res_summary$inter_p, method = "BH")
res_summary$g_ttest_pv_adj <- p.adjust(res_summary$g_ttest_pv, method = "BH")
nrow(res_summary[res_summary$inter_p_adj < 0.05,])
nrow(res_summary[res_summary$g_ttest_adj < 0.05,])
nrow(res_summary[res_summary$g_ttest_adj < 0.05 & res_summary$inter_p_adj < 0.05,])
paste0(row.names(res_summary[res_summary$g_ttest_adj > 0.05 & res_summary$inter_p_adj < 0.05,]), collapse = ",")

write.table(res_summary, file = paste0("summary_tables/", out_basepath, ".txt"), sep = "\t", quote = F, col.names = NA)

if (make_plots){
  dev.off()
}
