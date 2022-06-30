library(ggplot2)
library(RColorBrewer)
library('dplyr')
library('mgcv')
library(cluster)
library(factoextra)
library(psych)

source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/clustering_functions.R")


make_plots = F

n_points = 300
res_dif_all <- data.frame(age = seq(20, 75, length = n_points))
res_summary <- data.frame()

cnt = 1
if (make_plots){
  pdf(plot_path, width = 15, height = 15)
  par(mfrow=c(5,4)) 
}


indices = 1:ncol(traits_m)
#indices = sample(1:ncol(traits_m), 500, replace=F)
#plot_list <- list()

for (idx in indices){
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
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, T, 300, make_plots)

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

#####################################################################

#res_dif_all <- res_dif_all_genes[, c(2,rand)]

#res_dif_all <- read.delim("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.txt", header = T, sep = "\t", as.is = T, check.names = F)
#res_dif_all <- subset(res_dif_all, select = -1)
t_res_dif_all <- as.data.frame(t(res_dif_all))

grouped_res_dif_all <- data.frame(matrix(nrow = ncol(res_dif_all) - 1, ncol = 5))
colnames(grouped_res_dif_all) <- c("degree", "cluster", "indices", "intersection", "gene_names")
row.names(grouped_res_dif_all) <- colnames(res_dif_all)[-1]

for (g in colnames(res_dif_all)[-1]){
  #print(g)
  d <- get_poly_degree(res_dif_all, g)
  grouped_res_dif_all[g, "degree"] <- d
  
  int_pnts <- find_intersection_points(t_res_dif_all[g,], d)
  
  grouped_res_dif_all[g, "intersection"] <- paste(int_pnts, collapse = ",")
}

grouped_dif <- grouped_res_dif_all[order(grouped_res_dif_all$degree),]

grouped_dif$indices <- 0
gene_conv_table <- gene_table[match(row.names(grouped_dif), gene_table[,1]),]
all(gene_conv_table[,1] == row.names(grouped_dif))
grouped_dif$gene_names <- gene_conv_table[,2]
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
    num_k <- find_optimal_k(dif_gr, method = cl_method)
    clust <- do_clustering(dif_gr, num_k, method = cl_method, distance = cl_dist)
  }
  grouped_dif[grouped_dif$degree == d, "cluster"] = paste0(d, "_", clust)
}

grouped_dif$cluster <- paste(grouped_dif$cluster, grouped_dif$intersection, sep = "_")
grouped_dif <- grouped_dif[order(grouped_dif$cluster),]
#plot_clustering_ggplot(grouped_dif, paste0("../anova+gam/rnaseq_plots/ageDEgenes/new/LLD_expression_random_genes_clustering_intersection_", cl_method, "_", cl_dist, "_3.pdf"), plot_list)

plot_clustering(grouped_dif, paste0("../anova+gam/rnaseq_plots/ageDEgenes/new/LLD_expression_all_sign_genes_clustering_intersection_", cl_method, "_", cl_dist, ".pdf"))
write.table(grouped_dif, file = "../anova+gam/rnaseq_plots/ageDEgenes/new/LLD_expression_all_sign_genes_clustering_intersection_results.txt", sep = "\t", quote = F, col.names = NA)



#################
## Pathway enrichment
#################


library(enrichR)
dbs <- listEnrichrDbs()
db_sel <- c("BioPlanet_2019", "	KEGG_2019_Human", "Reactome_2016", "WikiPathways_2019_Human", "GO_Biological_Process_2018", "GO_Molecular_Function_2018")

sign_enrich <- data.frame(matrix(ncol = 5))
colnames(sign_enrich) <- c("cluster","db", "Term", "Overlap", "Adj_pval")
cnt <- 1
for (cl in unique(grouped_dif$cluster)){
  if (nrow(grouped_dif[grouped_dif$cluster == cl,]) > 5){
    enriched <- enrichr(grouped_dif[grouped_dif$cluster == cl, "gene_names"], db_sel)
    for (d in db_sel){
      r <- enriched[[d]][enriched[[d]]$Adjusted.P.value < 0.05,]
      if (nrow(r) > 0){
        for (rn in 1:nrow(r)){
          sign_enrich[cnt, "cluster"] <- cl
          sign_enrich[cnt, "db"] <- d
          sign_enrich[cnt, "Term"] <- r$Term[rn]
          sign_enrich[cnt, "Overlap"] <- r$Overlap[rn]
          sign_enrich[cnt, "Adj_pval"] <- r$Adjusted.P.value[rn]
          cnt = cnt + 1
        }
      }  
    }
  }
}


