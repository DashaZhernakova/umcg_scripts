
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/clustering_functions.R")


boxy = F

if (boxy){
  wd_path <- '/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/'
} else {
  wd_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/data/"
}
# Bile acids
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/BA37_50present_n1436.txt"
st_col = 1
# NMR metabolomics
traits_path <- "C:/Users/Dasha/work/UMCG/data/Metabolomics_shared_folder/2.metabolites/LLD_bloodlipids_nmr.txt"
st_col=5


gte_path <- "gte_all.txt"
pheno_path <- "age_gender_smk_contrac_cell_counts.cc_outliers_na.txt"
gene_table_path <- "geneid_to_gene_proteincoding_mainchr.txt"
plot_basepath <- "../anova+gam/nmr_lipidomics_plots.pdf"
correct_for_cellcounts = F
make_plots = T

setwd(wd_path)

# traits of interest
traits <- read.table(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- traits[,seq(st_col,ncol(traits))]

# Age, gender and other phenotypes
pheno0 <- as.data.frame(t(read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
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
res_dif_all <- data.frame(age = seq(20, 75, length = n_points))
res_summary <- data.frame()

cnt = 1
if (make_plots){
  pdf(plot_basepath, width = 15, height = 15)
  par(mfrow=c(5,4)) 
}


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
  if (length(trait_name) > 0){ #if gene id in gene convertion table
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
    
    res_dif = NULL
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, F, 300, make_plots, label = '')
    
    if (res_dif_lst[["inter_p"]] < 0.05){
      cnt <- cnt + 1
      res_dif_all[,trait_id] <- res_dif_lst[["dif"]]
      
      res_summary[trait_id,'inter_p'] = res_dif_lst[["inter_p"]]
      res_summary[trait_id,'g_beta'] = res_dif_lst[["g_beta"]]
      res_summary[trait_id,'g_pv'] = res_dif_lst[["g_pv"]]
      
      #res_summary[trait_id,'p1'] = res_dif_lst[["plots"]][[1]]
      #res_summary[trait_id,'p2'] = res_dif_lst[["p2"]]
      #plot_list[[idx]] = res_dif_lst[["plots"]]
    }
  }
  
}

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

grouped_dif[grouped_dif$degree == 2, "degree"] = 1 
grouped_dif[grouped_dif$degree > 3, "degree"] = 4 # set all degrees > 3 to 4, change later



cl_method = "pam"
cl_dist = "cor"
for (d in unique(grouped_dif$degree)){
  dif_gr <- t_res_dif_all[row.names(grouped_dif[grouped_dif$degree == d,]),]
  if (nrow(dif_gr) < 4){
    clust <- rep("1", nrow(dif_gr))
  } else {
    #num_k <- find_optimal_k(dif_gr, method = cl_method)
    num_k=3
    clust <- do_clustering(dif_gr, num_k, method = cl_method, distance = cl_dist)
  }
  grouped_dif[grouped_dif$degree == d, "cluster"] = paste0(d, "_", clust)
}

grouped_dif$cluster <- paste(grouped_dif$cluster, grouped_dif$intersection, sep = "_")
grouped_dif <- grouped_dif[order(grouped_dif$cluster),]
#plot_clustering_ggplot(grouped_dif, paste0("../anova+gam/rnaseq_plots/ageDEgenes/new/LLD_expression_random_genes_clustering_intersection_", cl_method, "_", cl_dist, "_3.pdf"), plot_list)

plot_clustering(grouped_dif, paste0("../anova+gam/NMR_clustering_k3_", cl_method, "_", cl_dist, ".pdf"))
#write.table(grouped_dif, file = "../anova+gam/LLD_expression_all_sign_genes_clustering_intersection_results.txt", sep = "\t", quote = F, col.names = NA)

