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
  #config_path <- "./config.yml"
  #script_folder <- "../../scripts/umcg_scripts/gender_differences/"
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
  pheno_table <- read.delim(paste0(config$basedir_path, "/", config$phenotype_table), sep = "\t", as.is = T, check.names = F)
  
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
add_inter_p_to_plot = config$add_inter_p_to_plot
plot_title = config$plot_title
outlier_correction = config$outlier_correction_method
outlier_correction_method <- config$outlier_correction_method
log_transform = config$log_transform
scale_transform = config$scale_transform
gam_family = config$gam_family
split_by_covariate = config$split_by_covariate
highlight_positive_in_split = config$highlight_positive_in_split
ttest_cutoff <- config$breakpoints_ttest_cutoff
deriv_cutoff <- config$breakpoints_derivates_cutoff
interp_cutoff <- ifelse("interp_cutoff" %in% names(config),  config$interp_cutoff, 0.05) 
write_fitted <- ifelse("write_fitted" %in% names(config),  config$write_fitted, F)
plot_points <- ifelse("plot_points" %in% names(config),  config$plot_points, T)
runCV <- ifelse("run_cross_validation" %in% names(config),  config$run_cross_validation, F)
ymax_hist <- ifelse("ymax_hist" %in% names(config),  config$ymax_hist, 1)
section_starts <- ifelse("sections" %in% names(config),  config$sections, character(0))

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

if (make_plots && config$plot_extention == "pdf"){
  if (nplotspp > 1){
    pdf(paste0(plot_path, ".pdf"), width = 15, height = 21)
    cat("Saving plots to ", plot_path, ".pdf", "\n")
     
  } else {
    pdf(paste0(plot_path, ".pdf"), width = 4, height = 5)
    cat("Saving plots to ", plot_path, ".pdf", "\n")
  }
} else if (make_plots && config$plot_extention == "png"){
  if (nplotspp > 1) {
    nrows <- ceiling(sqrt(num_traits))
    size <- 3*nrows
    png(paste0(plot_path, ".png"), width = size, height = size, units = 'in', res = 400)
    cat("Saving plots to ", plot_path, ".png", "\n")
    
  } else {
    png(paste0(plot_path, ".png"), width = 5, height = 4, units = 'in', res = 400)
    cat("Saving plots to ", plot_path, ".png", "\n")
  }
}


#
# Run the analyses
#
indices = 1:ncol(traits_m)
out_table_path <- paste0(out_basepath, "tables/", config$output_fname)

cnt = 1
cat("\nStarting the analyses\n")
for (idx in indices){
  if (config$plot_extention == "pdf" && cnt > nplotspp && make_plots){
    cnt = 1
    if (nplotspp > 1) par(mfrow=c(5,4))
  }
  trait_name <- ifelse(is.null(pheno_table), colnames(traits_m)[idx], pheno_table[pheno_table[,1] == colnames(traits_m)[idx], 2])
  cat(idx, " : ", trait_name, "\n")
  
  log_transform = FALSE
  if (colnames(traits_m)[idx] %in% pheno_to_log) log_transform = TRUE
  cat("\tLog tranform: ", log_transform, "\n")
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
  covariateslinear2 <- covariateslinear[covariateslinear != colnames(traits_m)[idx]]
  covariatesnonlinear2 <- covariatesnonlinear[covariatesnonlinear != colnames(traits_m)[idx]]
  
  res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, covariates_linear = covariateslinear2, covariates_nonlinear = covariatesnonlinear2, n_points = n_points, make_plots = make_plots, gam_family = gam_family, min_age = min_age, max_age = max_age, ymax_hist = ymax_hist, label = '', add_inter_p_to_plot = add_inter_p_to_plot, plot_title = plot_title, interp_cutoff = interp_cutoff, plot_points = plot_points, log_tr = log_transform, plot_density = T)
  if (res_dif_lst[["inter_p"]] < interp_cutoff){
    cnt <- cnt + 1
    cat("\tSignificant interaction detected.\n")
  }
    
}

  dev.off()


