library(cluster)
library(factoextra)
setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/")

script_folder <- "C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/"
source(paste0(script_folder, "/plotting_functions.R"))


fitted_matrix <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_summary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif_inters <- row.names(signif[signif$inter_p_adj_bonferroni < 0.05,])

d <- read.delim("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/results/tables/nmr_corrected_bmi_smk1_statins_bonferroni.gam_coefficients.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
#d <- read.delim("data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
#d <- na.omit(scale(d))
cormat <- cor(d, method = "spearman", use="pairwise.complete.obs")

num_k = 15
height = 1.5
#hc <- hclust(as.dist(1-cormat),  method = "average")
hc <- eclust(d, "hclust", hc_metric = "spearman", hc_method = "average")

#hc <- hclust(dist(d, method = 'euclidean'), method = "average")
#cl=cutree(hc, h = height)
cl = cutree(hc, k = num_k)
num_k <- length(unique(cl))
num_k

hc <- hkmeans(d,num_k,  hc.method = "average")
cl <- hc$cluster
clusters <- as.data.frame(cl)
clusters$signif_inter <- row.names(clusters) %in% signif_inters

ss <- silhouette(cl, as.dist(1-cormat))
mean_sil_score <- aggregate(ss[,3]~ss[,1], FUN=mean)

pdf(paste0("results/plots/nmr_clustering_hclust_", num_k, ".4.v2022.pdf"), width = 20, height = 20)
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
