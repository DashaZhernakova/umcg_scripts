d <- read.delim("data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)

fitted_matrix <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_all_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_summary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif_inters <- row.names(signif[signif$inter_p_adj_bonferroni < 0.05,])

d <- na.omit(d)
#d <- d[,colnames(d) %in% signif_inters]

#d <- d[,!colnames(d) %in% signif_inters]

num_k = 5
cormat <- cor(d, method = "spearman")
km <- kmeans(cor2dist(cormat), num_k, nstart = 50)
cl <- km$cluster

clusters <- as.data.frame(cl)
clusters$signif_inter <- row.names(clusters) %in% signif_inters
write.table(clusters, file = paste0("results/tables/nmr_clustering_kmeans.k",num_k,".txt"), sep = "\t", quote = F, col.names = NA)


pdf(paste0("../nmr_clustering_kmeans_", num_k, ".pdf"), width = 20, height = 20)
par(mfrow=c(3,3))
cnt <- 1
for (i in 1:max(cl)){
  lipids <- names(cl[cl == i])
  print(length(lipids))
  print(i)
  
  draw_multiple_fitted_lines(fitted_matrix[,lipids], signif_inters, plot_title = paste0("cluster ", cnt, "\nN = ", length(lipids)))
  cnt <- cnt + 1
}
dev.off()

