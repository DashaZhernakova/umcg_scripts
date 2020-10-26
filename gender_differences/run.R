
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/clustering_functions.R")
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/get_breakpoints.R")
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/calc_derivatives.R")

boxy = F

if (boxy){
  wd_path <- '/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/'
} else {
  wd_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/data/"
}
setwd(wd_path)

# Phenotypes
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/20170123_selection_phenotypes_for_TL_quant.txt"
st_col = 3

#Proteins
traits_path <- "CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt"
st_col = 1
traits0 <- as.data.frame(t(read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
out_basepath <- paste0("../plots_all_pheno/v2/proteins_with_f2")


# telomeres
traits_path <- "C:/Users/Dasha/work/UMCG/data/telomeres/telomeres_raw_t.txt"
st_col = 1
traits0 <- as.data.frame(t(read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

# Bile acids
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/BA37_50present_n1436.txt"
out_basepath <- paste0("../plots_all_pheno/v2/bile_acids_with_breakpoints_t3_d1.5e-4")

st_col = 1
# NMR metabolomics
traits_path <- "C:/Users/Dasha/work/UMCG/data/Metabolomics_shared_folder/2.metabolites/LLD_bloodlipids_nmr.txt"
st_col=5
out_basepath <- paste0("../plots_all_pheno/v2/NMR_with_f2")

# Untargeted metabolomics
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/data_1442samples_LLD_baseline_1183plasma_metabolites.txt"
st_col <- 1

#TMAO
traits_path <- "C:/Users/Dasha/work/UMCG/data/MR/results2/AA_T2D/data/raw_n1010/tmao.lld_raw_n1010.txt"
st_col <- 1
out_basepath <- paste0("../plots_all_pheno/v2/tmao_t4_d1.5e-4")

# Immune markers and cytokines
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/cytokines_and_immune_markers.txt"
st_col <- 7

#Smoking
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/Smoking_LLD_GoNL_for_Sasha.txt"
st_col <- 1
#gam_family = binomial(link = "logit")

# Diseases
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/LLD_diseases_1660.txt"
st_col <- 1
# gam_family = binomial(link = "logit")

# traits_m[traits_m[,"T2Dy_n"] == 2,] = 1
#if (length(unique(merged_tab[,1])) > 2){
#print("Skipping! Non binary")
#} else {

# Drugs
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/20150618_45drugs_1135metasubjects_atleast5users.txt"
st_col <- 1

# Diet
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/LLD_diet_1660.txt"
st_col <- 1
# merged_tab[,1] = merged_tab[,1] + 1
#gam_family = ocat(R = max(merged_tab[,1])

traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/diet_new_sept2015/Grams per item_group_ID.txt"
st_col = 1

traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/diet_new_sept2015/OV11_0098 Resultaten per deelnemer_ID.txt"
st_col = 2

traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/diet_new_sept2015/OV11_0098 Resultaten per item_ID.txt"

traits0 <- read.delim(traits_path, header = T, sep = "\t", as.is = T, check.names = F)

val_var <- 'ALCOHOL'
traits1 <- dcast(traits0, LLDEEPID_LL~ITEMNAAM, value.var = val_var)
row.names(traits1) <- traits1[,1]
traits <- sapply(traits1[,-1], function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits1)


###############################################################
gte_path <- "gte_all.txt"
pheno_path <- "age_gender_cell_counts_070920.txt"
gene_table_path <- "geneid_to_gene_proteincoding_mainchr.txt"

correct_for_cellcounts = F
make_plots = T
add_breakpoints = F

setwd(wd_path)

# traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits0 <- traits0[,seq(st_col,ncol(traits0))]
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)

#Age, gender and other phenotypes
pheno0 <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
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
ttest_cutoff <- 3
deriv_cutoff <- 0.00015
covariateslinear <- c("ba", "eo", "er", "gr", "ly", "mo", "tr")
covariatesnonlinear <- c()

out_basepath <- paste0("../plots_all_pheno/v2/cytokines_with_breakpoints_intervals_t3_d1.5e-4_v2")

res_dif_all <- data.frame(age = seq(min_age, max_age, length = n_points))
res_pdat_all <- data.frame(age = c(seq(min_age, max_age, length = n_points), seq(min_age, max_age, length = n_points)))
res_summary <- data.frame()
res_dif_lst <- data.frame()

cnt = 1
if (make_plots){
  pdf(paste0(out_basepath, ".pdf"), width = 15, height = 15)
  par(mfrow=c(5,4)) 
}


indices = 1:ncol(traits_m)
#indices = c(1,2,8,9,10,12,13,21,22,23,24,25) # cell type composition
cnt = 1

for (idx in indices){
  if (cnt > nplotspp & make_plots){
    cnt = 1
    par(mfrow=c(5,4))
  }
  print(idx)

  trait_id <- colnames(traits_m)[idx]
  trait_name = trait_id
  if (length(trait_name) > 0 & length(unique(traits_m[,idx])) > 1){ #if gene id in gene convertion table
    
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "zscore", log_tr = F, scale_tr = T)
    res_dif = NULL
    sex_dif_pval <- calculate_sex_diff_ttest(merged_tab, covariates = c(covariateslinear, covariatesnonlinear), min_age, max_age)
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, covariates_linear = covariateslinear, covariates_nonlinear = covariatesnonlinear, n_points = n_points, make_plots = make_plots, gam_family = gaussian(), label = '', add_breakpoints = add_breakpoints,  t_threshold = ttest_cutoff, derivatives_cutoff = deriv_cutoff)

    if (res_dif_lst[["inter_p"]] < 0.05){
      cnt <- cnt + 1
      res_dif_all[,trait_id] <- res_dif_lst[["dif"]]
      #res_pdat_all[,trait_id] <- res_dif_lst$pdat$pred
    }
      res_summary[trait_id,'inter_p'] = res_dif_lst[["inter_p"]]
      res_summary[trait_id,'g_beta'] = res_dif_lst[["g_beta"]]
      res_summary[trait_id,'g_pv'] = res_dif_lst[["g_pv"]]
      res_summary[trait_id,'g_ttest_pv'] = sex_dif_pval
      res_summary[trait_id,'cohen_f2'] = res_dif_lst[["cohen_f2"]]
      
      if (!is.null(res_dif_lst[['breakpoints']])){
        res_summary[trait_id, "breakpoints_men"] = ifelse(length(res_dif_lst[['breakpoints']][[2]]) > 0, res_dif_lst[['breakpoints']][[2]], "NA")
        res_summary[trait_id, "breakpoints_women"] = ifelse(length(res_dif_lst[['breakpoints']][[1]]) > 0, res_dif_lst[['breakpoints']][[1]], "NA")
      }
  }
}
res_summary$inter_p_adj <- p.adjust(res_summary$inter_p, method = "BH")
res_summary$g_ttest_pv_adj <- p.adjust(res_summary$g_ttest_pv, method = "BH")
write.table(res_summary, file = paste0(out_basepath, ".txt"), sep = "\t", quote = F, col.names = NA)
nrow(res_summary[res_summary$inter_p_adj < 0.05,])
nrow(res_summary[res_summary$g_ttest_pv_adj < 0.05,])
nrow(res_summary[res_summary$g_ttest_pv_adj < 0.05 & res_summary$inter_p_adj < 0.05,])
paste0(row.names(res_summary[res_summary$g_ttest_pv_adj > 0.05 & res_summary$inter_p_adj < 0.05,]), collapse = ",")
if (make_plots){
  dev.off()
}

###################################################
## Clustering
#

t_res_dif_all <- as.data.frame(t(res_dif_all))

grouped_res_dif_all <- data.frame(matrix(nrow = ncol(res_dif_all) - 1, ncol = 5))
colnames(grouped_res_dif_all) <- c("degree", "cluster", "indices", "intersection", "gene_names")
row.names(grouped_res_dif_all) <- colnames(res_dif_all)[-1]

for (g in colnames(res_dif_all)[-1]){
  #print(g)
  d <- get_poly_degree(res_dif_all, g)
  grouped_res_dif_all[g, "degree"] <- d
  
  #int_pnts <- find_intersection_points(t_res_dif_all[g,], d)
  #grouped_res_dif_all[g, "intersection"] <- paste(int_pnts, collapse = ",")
}

grouped_dif <- grouped_res_dif_all[order(grouped_res_dif_all$degree),]

grouped_dif$indices <- 0
grouped_dif$indices <- match(row.names(grouped_dif), colnames(traits_m), nomatch = 0)

#grouped_dif[grouped_dif$degree > 1, "degree"] = 2 
grouped_dif[grouped_dif$degree == 2, "degree"] = 1 
grouped_dif[grouped_dif$degree > 3, "degree"] = 4 # set all degrees > 3 to 4, change later



cl_method = "pam"
cl_dist = "cor"
for (d in unique(grouped_dif$degree)){
  dif_gr <- t_res_dif_all[row.names(grouped_dif[grouped_dif$degree == d,]),]
  if (nrow(dif_gr) < 4){
    clust <- rep("1", nrow(dif_gr))
  } else {
    num_k <- find_optimal_k(dif_gr, method = cl_method)
    #if (d == 1){
    #  num_k = 2
    #} else {
    #  num_k=4
    #}
    clust <- do_clustering(dif_gr, num_k, method = cl_method, distance = cl_dist)
  }
  grouped_dif[grouped_dif$degree == d, "cluster"] = paste0(d, "_", clust)
}

#grouped_dif$cluster <- paste(grouped_dif$cluster, grouped_dif$intersection, sep = "_")

# Cluster all together
num_k=6
clust <- do_clustering(t_res_dif_all[-1,], num_k, method = cl_method, distance = cl_dist)
grouped_dif[names(clust), "cluster"] <- clust


grouped_dif <- grouped_dif[order(grouped_dif$cluster),]
#plot_clustering_ggplot(grouped_dif, paste0("../anova+gam/rnaseq_plots/ageDEgenes/new/LLD_expression_random_genes_clustering_intersection_", cl_method, "_", cl_dist, "_3.pdf"), plot_list)

plot_clustering(grouped_dif, paste0("../anova+gam/metabolomics_clustering_", cl_method, "_", cl_dist, ".pdf"), scale_traits = T)
write.table(grouped_dif, file = paste0("../anova+gam/NMR_clustering_scaled_all_k6_", cl_method, "_", cl_dist, ".txt"), sep = "\t", quote = F, col.names = NA)


