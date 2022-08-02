library(reshape2)
library(dendextend)
library(rprojroot)
library(tidyverse)
library(circlize)

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

log_transform = config$log_transform
scale_transform = config$scale_transform

#
# start
#
indices = 1:ncol(traits_m)

res_rsq_table <- data.frame(matrix(nrow = length(indices), ncol = 4))
res_rsq_table_perc <- data.frame(matrix(nrow = length(indices), ncol = 4))
row.names(res_rsq_table) <- colnames(traits_m)[indices]
colnames(res_rsq_table) <- c("null_model", "sex", "age", "inter")
row.names(res_rsq_table_perc) <- colnames(traits_m)[indices]
colnames(res_rsq_table_perc) <- c("null_model", "sex", "age", "inter")
cnt = 1
cat("\nStarting the analyses\n")
for (idx in indices){
  
  trait_name <- ifelse(is.null(pheno_table), colnames(traits_m)[idx], pheno_table[pheno_table[,1] == colnames(traits_m)[idx], 2])
  cat(idx, " : ", trait_name, "\n")
  
  log_transform = FALSE
  if (colnames(traits_m)[idx] %in% pheno_to_log) log_transform = TRUE
  cat("\tLog tranform: ", log_transform, "\n")
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = log_transform, scale_tr = scale_transform)

  covariateslinear2 <- covariateslinear[covariateslinear != colnames(traits_m)[idx]]
  covariatesnonlinear2 <- covariatesnonlinear[covariatesnonlinear != colnames(traits_m)[idx]]
  
  exp_var_dif <- calc_explained_variance_CV(merged_tab, trait_name, covariates_linear = covariateslinear2, covariates_nonlinear = covariatesnonlinear2)
  exp_var_dif[exp_var_dif < 0] <- 0
  res_rsq_table[cnt,] <- exp_var_dif
  res_rsq_table_perc[cnt,] <- exp_var_dif/sum(exp_var_dif)
  cnt <- cnt + 1
}
write.table(res_rsq_table, file = paste0(out_basepath, "/tables/explained_variance_olink.txt"), sep = "\t", quote = F, row.names = F)
#res_rsq_table <- read.delim(paste0(out_basepath, "/tables/explained_variance_nmr.txt"), sep = "\t",  as.is = T, header = 1)
#row.names(res_rsq_table) <- colnames(traits_m)

 

# plot

max_h <- max(rowSums(res_rsq_table))
t <- (res_rsq_table/max_h)
t$rest <- 1-rowSums(t)
barcolor = c("#BDBDBD", "#66C2A5", "#FC8D62", "#8DA0CB",  "#FFFFFF")
labelcolor <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 12, name = "Paired"))

fitted_matrix <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
hc <- hclust(dist(t(fitted_matrix)), method = "average")

t["axis",] <- rep(0, ncol(t))

labels=c(hc$labels[hc$order], "axis")
ord <- c(hc$order, nrow(t))
dend <- as.dendrogram(hc)

circos.clear()
pdf(paste0(out_basepath, "plots/circos_barplot_nmr.gam_hclust_eucl_fitted_cv3.pdf"),width = 20, height = 20, useDingbats = F)
#png(paste0(out_basepath, "plots/circos_barplot_olink.abs.hclust_orig_cv.png"),width = 2000, height = 2000, res = 300)

circos.par(cell.padding=c(0,0,0,0), start.degree = 0)
circos.initialize("a",xlim=c(0,nrow(t)))

circos.track(ylim=c(0,1),track.height=.2,track.margin=c(0,0),bg.border=NA,
             panel.fun=function(x,y)for(i in 1:nrow(t))circos.text(i-.5,0,labels[i],adj=c(0,.5),
                                                                   facing="clockwise",niceFacing=T,cex=1.3,col= "black"))

circos.track(ylim=c(0,1),track.height=.3,track.margin=c(0,.01),bg.border=NA,
             panel.fun=function(x,y)circos.barplot(as.matrix(t)[ord,],-.5+1:nrow(t),
                                                   col=barcolor,bar_width=1,lwd=.3,border="lightgray"))

circos.track(ylim=c(0,attr(dend,"height")),track.height=.4,track.margin=c(0,0),
             bg.border=NA,panel.fun=function(x,y)circos.dendrogram(dend))

circos.clear()
dev.off()


#
# Plot selected trajectories
#
signif <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_summary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
ids <- pheno_table[match(row.names(signif), pheno_table[,2], nomatch = 0),1]
row.names(signif) <- ids

signif_inters <- row.names(signif[signif$inter_p_adj_bonferroni < 0.05,])

pdf(paste0("results/plots/nmr_trajectories_selected.v2022.pdf"), width = 20, height = 20)
par(mfrow=c(5,4))

all_clusters_fitted <- data.frame()
st <- "Crea"
end <- "Leu"
lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))

st <- "S.HDL.PL_p"
end <- "S.LDL.PL_p"

lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- t(fitted_matrix[,lipids])

st <- "XL.HDL.FC"
end <- "L.HDL.CE"

lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

st <- "S.HDL.C_p"
end <- "IDL.FC"
lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

st <- "XL.HDL.C"
end <- "TotCho"

lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

st <- "L.VLDL.TG"
end <- "M.HDL.TG_p"
lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

dev.off()

write.table(all_clusters_fitted, file =  "results/tables/nmr_trajectories_selected.v2022.txt", sep = "\t", quote = F, col.names = NA)

