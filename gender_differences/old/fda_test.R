library(funHDDC)
res_dif_all <- read.delim("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.txt", header = T, sep = "\t", as.is = T, check.names = F)
ages <- res_dif_all$age
res_dif_all <- subset(res_dif_all, select = -c(1,2))
row.names(res_dif_all) <- ages
t_res_dif_all <- as.matrix(t(res_dif_all))
nrow(t_res_dif_all)
ncol(t_res_dif_all)

num_k = 3
basis<- create.bspline.basis(c(0,1), nbasis=100)
var1<-smooth.basis(argvals=seq(0,1,length.out = 300),y=t(t_res_dif_all),fdParobj=basis)$fd
res.uni<-funHDDC(var1,K=3,model="AkBkQkDk",init="kmeans",threshold=0.2)
plot(var1,col=res.uni$class)

clusters <- data.frame(matrix(nrow = nrow(t_res_dif_all), ncol = 2))
colnames(clusters) <- c("cluster", "indices")
row.names(clusters) <- var1$fdnames["reps"][[1]]
clusters$cluster <- res.uni$class
clusters <- clusters[order(clusters$cluster),]

gene_conv_table <- gene_table[match(row.names(clusters), gene_table[,2]),]
all(gene_conv_table[,2] == row.names(clusters))
clusters$indices <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)

plot_path <- paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_FDA_",  num_k, "_2.diff.sorted.cor.pdf")
pdf(plot_path, width = 15, height = 15)
par(mfrow=c(5,4))

cur_clst <- 1
cnt <- 1

for (idx in clusters$indices){
  print (cnt)
  if (clusters[cnt,'cluster'] != cur_clst){
    cur_clst <- clusters[cnt,'cluster']
    
    par(mfrow=c(5,4))
  }
  if (idx != 0){
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, row.names(clusters)[cnt], correct_for_cellcounts, nrow(res_dif_all), T, label = paste0('cluster ', cur_clst))
  }
  cnt = cnt + 1
}
dev.off()


#
# Multivariate
#
gene_conv_table <- gene_table[match(row.names(t_res_dif_all), gene_table[,2]),]
ind <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)
res_table <- matrix(nrow = 2*300, ncol = length(ind))
colnames(res_table) <- gene_conv_table[,2]

for (idx in ind){
  trait_id <- colnames(traits_m)[idx]
  trait_name = gene_table[gene_table[,1] == trait_id,2]
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, correct_for_cellcounts, n_points, FALSE)
  res_table[,trait_name] <- res_dif_lst[["pdat"]]$pred

}
res_table = as.matrix(res_table)

num_k = 4
basis<- create.bspline.basis(c(0,1), nbasis=100)
var1<-smooth.basis(argvals=seq(0,1,length.out = 300),y=res_table[1:300,],fdParobj=basis)$fd
var2<-smooth.basis(argvals=seq(0,1,length.out = 300),y=res_table[301:nrow(res_table),],fdParobj=basis)$fd
res.multi<-funHDDC(list(var1, var2),K=3,model="AkBkQkDk",init="random",threshold=0.2)
plot(var1,col=res.uni$class)
