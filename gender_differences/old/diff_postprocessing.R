library(reshape2)
library(ggplot2)
library(pheatmap)
library(TSclust)
library(TSdist)
library(kmlShape)
library(psych)

res_dif_all <- read.delim("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.txt", header = T, sep = "\t", as.is = T, check.names = F)
ages <- res_dif_all$age
res_dif_all <- subset(res_dif_all, select = -c(1,2))
row.names(res_dif_all) <- ages
cor_matrix <- cor(res_dif_all)
#cor_matrix <- reorder_cormat(cor_matrix)
t_res_dif_all <- as.matrix(t(res_dif_all))

#
# TSclust & TSdist
#

#dist = "FRECHET"
#dist = "DTWARP"
dist = "PER"
pdf('../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.heatmap_COR.pdf', width = 15, height = 15, useDingbats = F)
#h <- pheatmap(cor_matrix, cluster_rows = T, cluster_cols = T, fontsize = 8)
D1 <- diss(t(res_dif_all), dist)

C1 <- hclust(D1)
h <- pheatmap(D1, cluster_rows = C1, cluster_cols = C1, labels_row = C1$labels[C1$order], labels_col = C1$labels[C1$order], clustering_method = 'complete', clustering_distance_rows = D1, clustering_distance_cols = D1, fontsize = 8)
dev.off()
clusters <- as.data.frame(sort(cutree(C1, k = num_k)))

num_k = 4

#
# kmlShape
#
#myClds <- cldsWide(t_res_dif_all, ages, row.names(t_res_dif_all))
#kmlShape(myClds, nbClusters = num_k)
#clusters <- as.data.frame(sort(myClds['clusters']))


colnames(clusters) <- c('cluster')
clusters$indices <- 0
gene_conv_table <- gene_table[match(row.names(clusters), gene_table[,2]),]
all(gene_conv_table[,2] == row.names(clusters))
clusters$indices <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)


plot_path <- paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_", dist, "_", num_k, ".diff.sorted.cor.pdf")
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


#https://datawookie.netlify.com/blog/2017/04/clustering-time-series-data/

 # reorder_cormat <- function(cormat){
 #   # Use correlation between variables as distance
 #   dd <- as.dist((1-cormat)/2)
 #   hc <- hclust(dd)
 #   cormat <-cormat[hc$order, hc$order]
 # }

cnt <- 1
n_pheno = length(clusters$indices)
step = 10
n_win = (70-20)/step + 1
diff_mw <- as.data.frame(matrix(nrow = n_win*n_pheno, ncol = 4))
colnames(diff_mw) <- c("pheno", "age", "men", "women")
for (idx in clusters$indices){ 
  print(cnt)
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  colnames(merged_tab)[1] <- "phenotype"
  women <- merged_tab$gender_F1M2 == 1
  avg_w <- mean(merged_tab[women, 1])
  avg_m <- mean(merged_tab[!women, 1])
  
  
  bin_num <- 1
  for (bin_st in seq(20, 70, step)){
    bin_end <- bin_st + step - 1
    
    bin_w <- merged_tab$age >= bin_st & merged_tab$age < bin_end & merged_tab$gender_F1M2 == 1
    bin_m <- merged_tab$age >= bin_st & merged_tab$age < bin_end & merged_tab$gender_F1M2 == 2
    diff_mw[cnt,1] <- pheno_name
    diff_mw[cnt,2] <- bin_st + step/2
    diff_mw[cnt,3] <- mean(merged_tab[bin_m, 1]) 
    diff_mw[cnt,4] <- mean(merged_tab[bin_w, 1])
    
    bin_num <- bin_num + 1
    cnt = cnt + 1
  }
  
}

ip <- getProfiles(t = "age", y = c("men", "women"), id = "pheno", data = diff_mw) 
plotProfiles(ip = ip, data = diff_mw, var = "men", tvar = "age")
mod <- GLMM_MCMC(y = diff_mw[,c("men", "women")], dist = c("gaussian", "gaussian"), 
                 id = diff_mw[,"pheno"], x = list(men = "empty", women = "empty"), z = list(men = diff_mw[, "age"], women = diff_mw[, "age"]),
                 random.intercept = c(TRUE, TRUE),
                 prior.b = list(Kmax = 2),
                 nMCMC = c(burn = 100, keep = 1000, thin = 10, info = 100),
                 parallel = FALSE)
mod <- NMixRelabel(mod, type = "stephens", keep.comp.prob = TRUE)




#
# kml3d
#
library("kml3d")

gene_conv_table <- gene_table[match(row.names(t_res_dif_all), gene_table[,2]),]
ind <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)
res_table <- matrix(ncol = 2*300, nrow = length(ind) + 1)
row.names(res_table) <- c("gender", gene_conv_table[,2])
res_table[1,] <- c(rep(1,300), rep(2, 300))
colnames(res_table) <- as.character(c(seq(20, 75, length = 300), seq(20, 75, length = 300)))
for (idx in ind){
  trait_id <- colnames(traits_m)[idx]
  trait_name = gene_table[gene_table[,1] == trait_id,2]
  print(idx)
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, correct_for_cellcounts, n_points, FALSE)
  
  pdat <- res_dif_lst[["pdat"]]
  pdat <- pdat[order(pdat$gender_F1M2, pdat$age),]
  
  res_table[trait_name,] <- as.numeric(pdat$pred)
  
}
res_table2 <- as.data.frame(res_table[1:30,])
#plot(x=as.numeric(colnames(res_table)), y=res_table["TFPI",], col = res_table["gender",])

(option1 <- parKml3d(scale = T))

cld3d_res <- cld3d(res_table2[-1,],timeInData=list(women=1:300,men=301:600), idAll = row.names(res_table2[-1,]))
kml3d(cld3d_res,3:4,nbRedrawing=20)

dist = 'kml3d'
num_k=3

clusters <- as.data.frame(matrix(ncol = 2, nrow = nrow(res_table2[-1,])))
clusters[,1] <- cld3d_res[paste0('c', num_k)][[1]]['clusters']
row.names(clusters) <- row.names(res_table2[-1,])
colnames(clusters) <- c('cluster', 'indices')

gene_conv_table <- gene_table[match(row.names(clusters), gene_table[,2]),]
all(gene_conv_table[,2] == row.names(clusters))
clusters$indices <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)





##########################################################
#
# Distances
#
x1 <- "ULK2" #increasing line
x2 <- "SAMHD1" #increasing line

x3 <- "STON2" #weaving line
x4 <- "TFPI" #weaving line

x5 <- "VSIG2" #half weaving line
x6 <- "LTBP1" #half weaving line

x7 <- "ASIC1" #decreasing line
x8 <- "FKBP1B" #decreasing line


t_res_dif_subs <- t_res_dif_all[c(x1,x2,x3,x4,x5,x6,x7,x8),]
dist_lst <- c('CDM', 'CID', 'COR', 'CORT', 'DTWARP', 'DWT','EUCL', 'INT.PER', 'NCD', 'PACF',  'PER',   'SPEC.LLR', 'SPEC.GLK', 'SPEC.ISD')
dist_matrix <- as.data.frame(matrix(ncol=10, nrow = length(dist_lst)))
row.names(dist_matrix) <- dist_lst

for (dist in dist_lst){
print(dist)
dist_matrix[dist, 1] <- as.matrix(diss(t_res_dif_subs, dist))[x1,x2]
dist_matrix[dist, 2] <- as.matrix(diss(t_res_dif_subs, dist))[x3, x4]
dist_matrix[dist, 3] <- as.matrix(diss(t_res_dif_subs, dist))[x5, x6]
dist_matrix[dist, 4] <- as.matrix(diss(t_res_dif_subs, dist))[x7, x8]

dist_matrix[dist, 5] <- as.matrix(diss(t_res_dif_subs, dist))[x2, x3]
dist_matrix[dist, 6] <- as.matrix(diss(t_res_dif_subs, dist))[x4, x5]
dist_matrix[dist, 7] <- as.matrix(diss(t_res_dif_subs, dist))[x6, x7]
dist_matrix[dist, 8] <- as.matrix(diss(t_res_dif_subs, dist))[x4, x7]
dist_matrix[dist, 9] <- as.matrix(diss(t_res_dif_subs, dist))[x2, x5]
}

dist = 'COR'
dist_matrix2[dist, 1] <- cor(as.matrix(t(t_res_dif_subs)))[x1,x2]
dist_matrix2[dist, 2] <- cor(as.matrix(t(t_res_dif_subs)))[x3, x4]
dist_matrix2[dist, 3] <- cor(as.matrix(t(t_res_dif_subs)))[x5, x6]
dist_matrix2[dist, 4] <- cor(as.matrix(t(t_res_dif_subs)))[x7, x8]

dist_matrix2[dist, 5] <- cor(as.matrix(t(t_res_dif_subs)))[x2, x3]
dist_matrix2[dist, 6] <- cor(as.matrix(t(t_res_dif_subs)))[x4, x5]
dist_matrix2[dist, 7] <- cor(as.matrix(t(t_res_dif_subs)))[x6, x7]
dist_matrix2[dist, 8] <- cor(as.matrix(t(t_res_dif_subs)))[x4, x7]
dist_matrix2[dist, 9] <- cor(as.matrix(t(t_res_dif_subs)))[x2, x5]



####
##########
####

calcDist_plot <- function(d, i, j, make_plots){

  
  dist_val <- distFrechet(as.numeric(colnames(d)), d[i,], as.numeric(colnames(d)), d[j,])
  dist_val2 <- cor(d[i,], d[j,], method = "spearman")
  
  if (make_plots){
    ymax = max(d[i,], d[j,])
    ymin = min(d[i,], d[j,])
    plot(x=as.numeric(colnames(t_res_dif_subs)), y=d[i,], type='l', col = 'red', 
       ylim = c(ymin, ymax), main = paste(row.names(d)[i], "-", row.names(d)[j], "\n", dist_val, "\n", dist_val2))
    lines(x=as.numeric(colnames(d)), y=d[j,],  type='l', col = 'blue')
  }
  return(dist_val2)
}



dist_tab <- matrix(ncol = 3, nrow = nrow(t_res_dif_subs)*nrow(t_res_dif_subs))
cnt = 1
visited <- c()
for (i in 1:nrow(t_res_dif_subs)){
  visited <- c(visited, i)
  for (j in 1:nrow(t_res_dif_subs)){
    if (! j %in% visited){
    dist_tab[cnt, 1] = i
    dist_tab[cnt, 2] = j
    dist_tab[cnt, 3] = calcDist_plot(t_res_dif_subs, i, j, FALSE)
    cnt = cnt + 1
    }
  }
}
dist_tab <- dist_tab[order(dist_tab[,3]),]

pdf("../anova+gam/rnaseq_plots/ageDEgenes/test_dist_Frechet_cor_spearman2.pdf", width = 15, height = 15)
par(mfrow=c(5,5))
visited <- c()
for (r in 1:nrow(dist_tab)){
  i = dist_tab[r,1]
  j = dist_tab[r,2]

  g1 <- row.names(t_res_dif_subs)[i]
  g2 <- row.names(t_res_dif_subs)[j]
  g_id1 <- gene_table[gene_table[,2] == g1,1]
  idx1 <- which(colnames(traits_m) == g_id1)
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx1)
  res <- plot_scatter_and_gam2(merged_tab, g1, correct_for_cellcounts, n_points, T, label = "")
  g_id2 <- gene_table[gene_table[,2] == g2,1]
  idx2 <- which(colnames(traits_m) == g_id2)
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx2)
  res <- plot_scatter_and_gam2(merged_tab, g2, correct_for_cellcounts, n_points, T, label = "")
  calcDist_plot(t_res_dif_subs, i, j, TRUE)

}
dev.off()



#
###################################
#
# DBSCAN
##
library(dbscan)
cor_dist <- cor2dist(cor(res_dif_all))
dbscan::kNNdistplot(cor_dist, k =  10)
db<-dbscan(cor_dist,eps=0.8)
db
clusters <- as.data.frame(matrix(ncol = 2, nrow = nrow(cor_dist)))
clusters[,1] <- db$cluster
row.names(clusters) <- row.names(cor_dist)
colnames(clusters) <- c('cluster', 'indices')

gene_conv_table <- gene_table[match(row.names(clusters), gene_table[,2]),]
all(gene_conv_table[,2] == row.names(clusters))
clusters$indices <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)
clusters <- clusters[order(clusters$cluster),]
dist='dbscan'
num_k = 0.8


plot_path <- paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_", dist, "_", num_k, ".diff.sorted.cor.pdf")
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
