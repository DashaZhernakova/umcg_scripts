library(rprojroot)
library(tidyverse)
setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/")

script_folder <- "C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/"
config_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/config_local.yml"

cat("script folder:", script_folder, "\n")
source(paste0(script_folder, "/preprocessing_gam_fitting_functions.R"))
source(paste0(script_folder, "/get_breakpoints.R"))
source(paste0(script_folder, "/additional_functions.R"))
source(paste0(script_folder, "/plotting_functions.R"))
source(paste0(script_folder, "/clustering_functions.R"))

cat("Using config file: ", config_path, "\n")
config <- config::get(file = config_path)
# save the config in results folder
file.copy(config_path, paste0(config$basedir_path, "configs/", config$output_fname, "_cfg.yml"), overwrite = T)

#
# Read data
#

traits_path <- paste0(config$basedir_path, "/", config$traits_path)
pheno_path <- paste0(config$basedir_path, "/", config$pheno_path)
cat("Data paths:\nphenotype traits:", traits_path, "\r\ncovariates:", pheno_path, "\noutput base folder:", config$basedir_path, "\n\n")

# read phenotype traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)
traits2use <- unlist(strsplit(config$traits2use, ",")) # choose phenotypes to run the analysis for
if (length(traits2use) > 0) {
  traits <- as.data.frame(traits[,traits2use, drop = F])
  cat("Running the analysis only for a subset of phenotypes: ", paste(traits2use, collapse = ", "), "\n")
}


# read age, gender and other covariate phenotypes
pheno0 <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)


# Covariates
covariateslinear <- unlist(strsplit(config$covariateslinear, ","))
covariatesnonlinear <- unlist(strsplit(config$covariatesnonlinear, ","))

if (length(covariateslinear) > 0) print(paste0("covariates to add as linear terms in the gam model:", paste(covariateslinear, collapse = ", ")))
if (length(covariatesnonlinear) > 0) print(paste0("covariates to add as spline non-linear terms in the gam model:", paste(covariatesnonlinear, collapse = ", ")))


phenos2use <- unlist(strsplit(config$phenos2use, ","))
if (length(phenos2use) > 0) {
  pheno0 <- pheno0[,c("age", "gender_F1M2", phenos2use)] #choose covariate phenotypes to select from the file
} else {
  pheno0 <- pheno0[,c("age", "gender_F1M2", covariateslinear, covariatesnonlinear)]
}
pheno <- na.omit(pheno0)

#order samples in the two tables
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)

#traits_m <- traits_m[order(pheno_m$age), , drop = F]
#pheno_m <- pheno_m[order(pheno_m$age),]

cat("Number of available phenotypes: ", num_traits, "\n")
cat("Number of shared samples: ", nrow(traits_m), "\n")

covariates_before <- unlist(strsplit(config$covariates_before, ","))
if (length(covariates_before) > 0){
  print(paste0("Correcting for covariates using linear regression before gam fitting: ", paste(covariates_before, collapse = ", ")))
  traits_m <- correct_for_covariates_before(traits_m, pheno_m, covariates_before)
}

pheno_table <- NULL
if ("phenotype_table" %in% names(config)){
  pheno_table <- read.delim(config$phenotype_table, sep = "\t", as.is = T, check.names = F)
  
}

#
# Other parameters
#

nplotspp = config$n_plots_ppage
n_points = config$n_age_points
min_age = config$min_age
max_age = config$max_age
make_plots = config$make_plots
add_inter_p_to_plot = config$add_inter_p_to_plot
plot_title = config$plot_title
outlier_correction = config$outlier_correction_method
outlier_correction_method <- config$outlier_correction_method
log_transform = config$log_transform
scale_transform = config$scale_transform
interp_cutoff <- ifelse("interp_cutoff" %in% names(config),  config$interp_cutoff, 0.05) 
write_fitted <- ifelse("write_fitted" %in% names(config),  config$write_fitted, F)
plot_points <- ifelse("plot_points" %in% names(config),  config$plot_points, T)
runCV <- ifelse("run_cross_validation" %in% names(config),  config$run_cross_validation, F)
ymax_hist <- ifelse("ymax_hist" %in% names(config),  config$ymax_hist, 1)

if ("pheno_to_log" %in% names(config)){
  pheno_to_log <- unlist(strsplit(config$pheno_to_log, ","))
} else {
  pheno_to_log <- character(0)
}
cat("Phenotypes to log-transform: ", pheno_to_log, "\n")
#
# Plot initialization
#

out_basepath <- paste0(config$basedir_path, "/results/")
plot_path <- paste0(out_basepath, "plots/", config$output_fname)


#
# Run the analyses
#
indices = 1:ncol(traits_m)
out_table_path <- paste0(out_basepath, "tables/", config$output_fname)


age_breaks = c(20, 40, 60, 80)

lm_estimates <- data.frame(matrix(nrow = num_traits, ncol =  2*(length(age_breaks) -1)))
#lm_estimates <- data.frame(matrix(nrow = num_traits, ncol =  40))

cnt = 1
cat("\nStarting the analyses\n")
for (idx in indices){
  
  trait_name <- colnames(traits_m)[idx]
  cat(idx, " : ", trait_name, "\n")
  
  log_transform = FALSE
  if (colnames(traits_m)[idx] %in% pheno_to_log) log_transform = TRUE
  cat("\tLog tranform: ", log_transform, "\n")
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "No", log_tr = log_transform, scale_tr = scale_transform)
  lm_estimates[idx,] <- get_lm_estimates_v3(merged_tab, age_breaks = age_breaks)
  #lm_estimates[idx,] <- get_gam_summary(merged_tab, trait_name, covariateslinear, covariatesnonlinear )
  rownames(lm_estimates)[idx] <- trait_name
}
write.table(lm_estimates, paste0(out_table_path, ".gam_coefficients.txt"), sep = "\t", quote = F, col.names = NA)
# 
 nmr_corr <- cor(traits_m, use="complete.obs")
 corr <- cor(t(lm_estimates))
 test = hclust(dist(corr))
 ORDER = rownames(corr)[test$order]
 nmr_corr[ORDER,ORDER] ->  nmr_corr
 corr[ORDER,ORDER] ->  corr
 New_matrix = nmr_corr
 New_matrix[lower.tri(nmr_corr)] <- as.matrix(corr)[lower.tri(as.matrix(corr))]
# 
 pdf("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/results/plots/NMR_corrplot_test2.pdf", height = 40, width = 40)
 corrplot(New_matrix, diag=FALSE, tl.col="black", tl.cex = 0.5, method = 'color') -> Plot_m
 corrMatOrder(corr = Plot_m[test$order], order = names(corr)[test$order]  )
 dev.off()
# 
# pdf("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/results/plots/NMR_clustering.pdf", width = 15, height = 21)
# par(mfrow=c(5,4)) 
# 
# for (trait_name in ORDER){
#   if (config$plot_extention == "pdf" && cnt > nplotspp && make_plots){
#     cnt = 1
#     if (nplotspp > 1) par(mfrow=c(5,4))
#   }
#   idx = which (colnames(traits_m) == trait_name)
#   cat(trait_name, cnt, "\n", sep = " ")
#   merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
#   res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name,  make_plots = T , plot_title = trait_name)
#   cnt <- cnt + 1
# }
# dev.off()





fitted_matrix <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_all_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_summary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif_inters <- row.names(signif[signif$inter_p_adj_bonferroni < 0.05,])
num_k = 9
km <- kmeans(lm_estimates, num_k, nstart = 50)
cl <- km$cluster

clusters <- as.data.frame(cl)
clusters$signif_inter <- row.names(clusters) %in% signif_inters

ss <- silhouette(km$cluster, dist(cor(t(lm_estimates))))
mean_sil_score <- aggregate(ss[,3]~ss[,1], FUN=mean)

pdf(paste0("results/plots/nmr_clustering_hclust_", num_k, ".2.v2022.pdf"), width = 20, height = 20)
par(mfrow=c(5,4))
cnt <- 1
for (i in 1:max(cl)){
  lipids <- row.names(clusters[clusters$cl == i,])
  print(length(lipids))
  print(i)
  
  draw_multiple_fitted_lines(fitted_matrix[,lipids], signif_inters, plot_title = paste0("cluster ", cnt, "\nN = ", length(lipids)))
  cnt <- cnt + 1
}
dev.off()


for (i in c(1,3,2,5,6,4)){
  lipids <- row.names(clusters[clusters$cl == i,])
  print(length(lipids))
  print(i)
  
  draw_multiple_fitted_lines(fitted_matrix[,lipids], signif_inters, plot_title = paste0("cluster ", cnt, "\nN = ", length(lipids)))
  cnt <- cnt + 1
}
