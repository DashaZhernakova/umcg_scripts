library(rprojroot)
library(tidyverse)


getCurrentFileLocation <-  function()
{
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col=value, into=c("key", "value"), sep="=", fill='right') %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  if (length(this_file)==0)
  {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}

isRStudio <- Sys.getenv("RSTUDIO") == "1"
if (isRStudio) {
  script_folder <- "C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/"
  config_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/config_local.yml"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- args[1]
  script_folder <- getCurrentFileLocation()
}
cat("script folder:", script_folder, "\n")
source(paste0(script_folder, "/preprocessing_gam_fitting_functions.R"))
source(paste0(script_folder, "/get_breakpoints.R"))
source(paste0(script_folder, "/additional_functions.R"))
source(paste0(script_folder, "/plotting_functions.R"))

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
  traits <- traits[,traits2use]
  cat("Running the analysis only for a subset of phenotypes: ", paste(traits2use, collapse = ", "), "\n")
}


# read age, gender and other covariate phenotypes
pheno0 <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
phenos2use <- unlist(strsplit(config$phenos2use, ","))
if (length(phenos2use) > 0) pheno0 <- pheno0[,c("age", "gender_F1M2", phenos2use)] #choose covariate phenotypes to select from the file
pheno <- na.omit(pheno0)

#order samples in the two tables
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

cat("Number of available phenotypes: ", num_traits, "\n")
cat("Number of shared samples: ", nrow(traits_m), "\n")

# Covariates
covariateslinear <- unlist(strsplit(config$covariateslinear, ","))
covariatesnonlinear <- unlist(strsplit(config$covariatesnonlinear, ","))

if (length(covariateslinear) > 0) print(paste0("covariates to add as linear terms in the gam model:", paste(covariateslinear, collapse = ", ")))
if (length(covariatesnonlinear) > 0) print(paste0("covariates to add as spline non-linear terms in the gam model:", paste(covariateslinear, collapse = ", ")))

covariates_before <- unlist(strsplit(config$covariates_before, ","))
if (length(covariates_before) > 0){
  print(paste0("Correcting for covariates using linear regression before gam fitting: ", paste(covariates_before, collapse = ", ")))
  traits_m <- correct_for_covariates_before(traits_m, pheno_m, covariates_before)
}

#
# Other parameters
#

nplotspp = config$n_plots_ppage
n_points = config$n_age_points
min_age = config$min_age
max_age = config$max_age
make_plots = config$make_plots
add_breakpoints = config$add_breakpoints
outlier_correction = config$outlier_correction_method
outlier_correction_method <- config$outlier_correction_method
log_transform = config$log_transform
scale_transform = config$scale_transform
gam_family = config$gam_family
split_by_covariate = config$split_by_covariate
ttest_cutoff <- config$breakpoints_ttest_cutoff
deriv_cutoff <- config$breakpoints_derivates_cutoff

#
# Plot initialization
#

out_basepath <- paste0(config$basedir_path, "/results/")
plot_path <- paste0(out_basepath, "plots/", config$output_fname)

if (make_plots && config$plot_extention == "pdf"){
  if (nplotspp > 1){
    pdf(paste0(plot_path, ".pdf"), width = 15, height = 15)
    cat("Saving plots to ", plot_path, ".pdf", "\n")
    par(mfrow=c(5,4)) 
  } else {
    pdf(paste0(plot_path, ".pdf"), width = 5, height = 4)
    cat("Saving plots to ", plot_path, ".pdf", "\n")
  }
} else if (make_plots && config$plot_extention == "png"){
  if (nplotspp > 1) {
    nrows <- ceiling(sqrt(num_traits))
    size <- 3*nrows
    png(paste0(plot_path, ".png"), width = size, height = size, units = 'in', res = 400)
    cat("Saving plots to ", plot_path, ".png", "\n")
    par(mfrow=c(nrows,nrows))
  } else {
    png(paste0(plot_path, ".png"), width = 5, height = 4, units = 'in', res = 400)
    cat("Saving plots to ", plot_path, ".png", "\n")
  }
}


#
# Run the analyses
#

res_summary <- data.frame()
res_dif_lst <- data.frame()
out_table_path <- paste0(out_basepath, "tables/", config$output_fname)
cat("\nStarting the analyses\n")
indices = 1:ncol(traits_m)
cnt = 1

for (idx in indices){
  if (config$plot_extention == "pdf" && cnt > nplotspp && make_plots){
    cnt = 1
    if (nplotspp > 1) par(mfrow=c(5,4))
  }
  trait_name <- colnames(traits_m)[idx]
  cat(idx, " : ", trait_name, "\n")
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
  if (split_by_covariate == ""){
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, covariates_linear = covariateslinear, covariates_nonlinear = covariatesnonlinear, n_points = n_points, make_plots = make_plots, gam_family = gam_family, label = '', add_breakpoints = add_breakpoints,  t_threshold = ttest_cutoff, derivatives_cutoff = deriv_cutoff)
  } else {
    run_for_split_by_covariate(merged_tab, trait_name, covariate_to_split = split_by_covariate , covariates_linear = covariateslinear, covariates_nonlinear = covariatesnonlinear, n_points = n_points, make_plots = make_plots, gam_family = gam_family)
  }
  if (res_dif_lst[["inter_p"]] < 0.05){
    cnt <- cnt + 1
    cat("\tSignificant interaction detected.\n")
  }
  
  sex_dif_pval <- calculate_sex_diff_ttest(merged_tab, covariates = c(covariateslinear, covariatesnonlinear), min_age, max_age)
  res_summary[trait_name,'inter_p'] = res_dif_lst[["inter_p"]]
  res_summary[trait_name,'g_beta'] = res_dif_lst[["g_beta"]]
  res_summary[trait_name,'g_pv'] = res_dif_lst[["g_pv"]]
  res_summary[trait_name,'g_ttest_pv'] = sex_dif_pval
  res_summary[trait_name,'cohen_f2'] = res_dif_lst[["cohen_f2"]]
  
  if (!is.null(res_dif_lst[['breakpoints']])){
    res_summary[trait_name, "breakpoints_men"] = ifelse(length(res_dif_lst[['breakpoints']][[2]]) > 0, res_dif_lst[['breakpoints']][[2]], "NA")
    res_summary[trait_name, "breakpoints_women"] = ifelse(length(res_dif_lst[['breakpoints']][[1]]) > 0, res_dif_lst[['breakpoints']][[1]], "NA")
  }
}
res_summary$inter_p_adj <- p.adjust(res_summary$inter_p, method = "BH")
res_summary$g_ttest_pv_adj <- p.adjust(res_summary$g_ttest_pv, method = "BH")
write.table(res_summary, file = paste0(out_table_path, "_summary.txt"), sep = "\t", quote = F, col.names = NA)

if (make_plots){
  dev.off()
}

