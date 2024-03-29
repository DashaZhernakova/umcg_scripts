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
  script_folder <- "/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/"
  config_path <- "/Users/Dasha/work/UMCG/data/gender_differences/omics/results/config_local.yml"
} else {
  args <- commandArgs(trailingOnly = TRUE)
  config_path <- args[1]
  script_folder <- getCurrentFileLocation()
  #config_path <- "./config.yml"
  #script_folder <- "../../scripts/umcg_scripts/gender_differences/"
}
cat("script folder:", script_folder, "\n")
source(paste0(script_folder, "/preprocessing_gam_fitting_functions.R"))
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
plot_density <- ifelse("plot_density" %in% names(config), config$plot_density, FALSE)

if ("pheno_to_log" %in% names(config)){
   pheno_to_log <- unlist(strsplit(config$pheno_to_log, ","))
} else {
   pheno_to_log <- character(0)
}
cat("Phenotypes to log-transform: ", pheno_to_log, "\n")
#
# Plot initialization
#
if (runCV) make_plots = F

out_basepath <- paste0(config$basedir_path, "/results/")
plot_path <- paste0(out_basepath, "plots/", config$output_fname)

if (make_plots && config$plot_extention == "pdf"){
  if (nplotspp > 1){
    pdf(paste0(plot_path, ".pdf"), width = 15, height = 21)
    cat("Saving plots to ", plot_path, ".pdf", "\n")
    par(mfrow=c(5,4)) 
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
    par(mfrow=c(nrows,nrows))
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


if (runCV){
  cat("\nRunning only cross-validation\n")
  res_summary <- data.frame(matrix(nrow = num_traits, ncol = 4))
  colnames(res_summary) <- c("gam_inter", "gam_nointer", "lm_inter", "lm_nointer")
  rownames(res_summary) <- colnames(traits_m)
  for (idx in indices){
    trait_name <- ifelse(is.null(pheno_table), colnames(traits_m)[idx], pheno_table[pheno_table[,1] == colnames(traits_m)[idx], 2])
    cat(idx, " : ", colnames(traits_m)[idx], " : ", trait_name, "\n")
    log_transform = FALSE
    if (colnames(traits_m)[idx] %in% pheno_to_log) log_transform = TRUE
    cat("\tLog tranform: ", log_transform, "\n")
    
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
    cv_lst <- run_cross_validation(merged_tab, pheno_name, covariateslinear, covariatesnonlinear, min_age, max_age)
    #print (unlist(cv_lst))
    res_summary[trait_name,] <- unlist(cv_lst)
  }
  write.table(res_summary, file = paste0(out_table_path, "_cross_validation.txt"), sep = "\t", quote = F, col.names = NA)
  quit()
}

res_summary <- data.frame()
res_dif_lst <- data.frame()
if (write_fitted){
  fitted_lines <- data.frame(matrix(nrow = n_points*2, ncol = length(indices)))
  if (split_by_covariate == ""){
    fitted_res_table <- data.frame(matrix(ncol = 6, nrow = 600*num_traits))
    colnames(fitted_res_table) <- c("Phenotype", "age", "sex_F1M2", "predicted_value", "lower_CI", "upper_CI")
  } else {
    fitted_res_table <- data.frame(matrix(ncol = 7, nrow = 2*600*num_traits))
    colnames(fitted_res_table) <- c("Phenotype", "sex_F1M2", "covariate", "age", "predicted_value", "lower_CI", "upper_CI")
  }
}

cnt = 1
cnt_sign <- 1
plots <- list()

cat("\nStarting the analyses\n")
for (idx in indices){
  if (config$plot_extention == "pdf" && cnt > nplotspp && make_plots && plot_density == F){
    cnt = 1
    if (nplotspp > 1) par(mfrow=c(5,4))
  }
  trait_name <- ifelse(is.null(pheno_table), colnames(traits_m)[idx], pheno_table[pheno_table[,1] == colnames(traits_m)[idx], 2])
  cat(idx, " : ", trait_name, "\n")
  
  log_transform = FALSE
  if (colnames(traits_m)[idx] %in% pheno_to_log) log_transform = TRUE
  cat("\tLog tranform: ", log_transform, "\n")
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
  if (split_by_covariate == ""){
    new_covars <- remove_related_covariates(traits_m, covariateslinear, covariatesnonlinear)
    
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, covariates_linear = new_covars$linear, covariates_nonlinear = new_covars$nonlinear, n_points = n_points, make_plots = make_plots, gam_family = gam_family, min_age = min_age, max_age = max_age, ymax_hist = ymax_hist, label = '', add_inter_p_to_plot = add_inter_p_to_plot, plot_title = plot_title, interp_cutoff = interp_cutoff, plot_points = plot_points, add_breakpoints = add_breakpoints,  t_threshold = ttest_cutoff, derivatives_cutoff = deriv_cutoff, log_tr = log_transform, plot_density = plot_density)
    
    if (res_dif_lst[["inter_p"]] < interp_cutoff){
      if (plot_density && ! is.null(res_dif_lst$p)) plots[[cnt_sign]] <- res_dif_lst$p
      cnt <- cnt + 1
      cnt_sign <- cnt_sign + 1
      cat("\tSignificant interaction detected.\n")
    }
    
    res_summary[trait_name,'inter_p'] = res_dif_lst[["inter_p"]]
    res_summary[trait_name,'g_beta'] = res_dif_lst[["g_beta"]]
    res_summary[trait_name,'g_pv'] = res_dif_lst[["g_pv"]]
    res_summary[trait_name,'cohen_f2'] = res_dif_lst[["cohen_f2"]]
    
    sex_dif_lm <- calculate_sex_diff_lm(merged_tab, covariates = c(new_covars$linear, new_covars$nonlinear), min_age, max_age)
    res_summary[trait_name,'g_lm_pv'] = sex_dif_lm[[1]]
    res_summary[trait_name,'g_lm_f2'] = sex_dif_lm[[2]]
    
    age_sex_only <- fit_gam_age_gender_only(merged_tab, trait_name, covariates_linear = new_covars$linear, covariates_nonlinear = new_covars$nonlinear, n_points = n_points, gam_family = gam_family, min_age = min_age, max_age = max_age, log_tr = log_transform)
    res_summary[trait_name,'sex_gam_pv'] = age_sex_only[1]
    res_summary[trait_name,'sex_gam_f2'] = age_sex_only[2]
    res_summary[trait_name,'age_gam_pv'] = age_sex_only[3]
    res_summary[trait_name,'age_gam_f2'] = age_sex_only[4]

    
    if (write_fitted) {
      fitted_lines[,trait_name] = res_dif_lst[["pdat"]]$pred
      fitted_res_table <- write_fitted_data (res_dif_lst[["pdat"]], fitted_res_table, idx, trait_name)
    }

  } else {
    pdat_lst <- run_for_split_by_covariate(merged_tab, trait_name, covariate_to_split = split_by_covariate , highlight_positive = highlight_positive_in_split, covariates_linear = covariateslinear, covariates_nonlinear = covariatesnonlinear, n_points = n_points, make_plots = make_plots, gam_family = gam_family, plot_points = plot_points, log_tr = log_transform, min_age = min_age, max_age = max_age)
    if (write_fitted) fitted_res_table <- write_fitted_data (pdat_lst, fitted_res_table, idx, trait_name, split_by_covariate)
    #cat(cor(res_diff[seq(1,300)], res_diff[seq(301,600)], method = "spearman"))
  }
}
res_summary$inter_p_adj_BH <- p.adjust(res_summary$inter_p, method = "BH")
res_summary$g_lm_pv_adj_BH <- p.adjust(res_summary$g_lm_pv, method = "BH")
res_summary$age_gam_pv_adj_BH <- p.adjust(res_summary$age_gam_pv, method = "BH")
res_summary$sex_gam_pv_adj_BH <- p.adjust(res_summary$sex_gam_pv, method = "BH")


res_summary$inter_p_adj_bonferroni <- p.adjust(res_summary$inter_p, method = "bonferroni")
res_summary$g_lm_pv_adj_bonferroni <- p.adjust(res_summary$g_lm_pv, method = "bonferroni")
write.table(res_summary, file = paste0(out_table_path, "_summary.txt"), sep = "\t", quote = F, col.names = NA)
if (write_fitted) {
  write.table(fitted_lines, file = paste0(out_table_path, "_fitted.txt"), sep = "\t", quote = F, col.names = NA)
  write.table(fitted_res_table, file = paste0(out_table_path, "_datasource.txt"), sep = "\t", quote = F, col.names = NA)
}
if (plot_density) write_plots_cowplot(plots, plot_path)

if (make_plots){
  dev.off()
}
