

get_poly_degree <- function (g){
  x <- res_dif_all[,c("age",g)]
  lm0 <- chooseBestModel(x)
  
  if ("I(age^5)" %in% names(lm0$coefficients)){
    poly_d = 5
  } else if ("I(age^4)" %in% names(lm0$coefficients)){
    poly_d = 4
  } else if ("I(age^3)" %in% names(lm0$coefficients)) {
    poly_d = 3
  } else if ("I(age^2)" %in% names(lm0$coefficients)) {
    poly_d = 2
  } else if ("age" %in% names(lm0$coefficients)) {
    poly_d = 1
  } else {
    poly_d = 0
  }
  return (poly_d)
}
pdat_lm <- with(x, expand.grid(age = seq(20, 75, length = 500)))
pdat_lm <- transform(pdat_lm, pred = predict(lm0, newdata = pdat_lm, type = "response"))

plot(x)
lines(pred ~ age, data = pdat_lm, col = "red", lwd = 2)

chooseBestModel <- function(d){
  colnames(d)[2] = "pheno"
  lm1 <- lm(pheno~age, data = d)
  if (summary(lm1)$adj.r.squared > 0.98) return (lm1)
  lm2 <- lm(pheno~age + I(age^2) , data = d)
  if (summary(lm2)$adj.r.squared > 0.98) return (lm2)
  lm3 <- lm(pheno~age + I(age^2) + I(age^3) , data = d)
  if (summary(lm3)$adj.r.squared > 0.98) return (lm3)
  lm4 <- lm(pheno~age + I(age^2) + I(age^3) + I(age^4), data = d)
  if (summary(lm4)$adj.r.squared > 0.98) return (lm4)
  lm5 <- lm(pheno~age + I(age^2) + I(age^3) + I(age^4) + I(age^5), data = d)
  if (summary(lm5)$adj.r.squared > 0.98) return (lm5)
  
  models <- list(lm1, lm2, lm3, lm4, lm5)
  AICs = NULL
  BICs = NULL
  for(i in 1:length(models)){
    AICs[i] = AIC(models[[i]])
    BICs[i] = BIC(models[[i]])
  }
  
  #see which models were chosen as best by both methods
  aic_n <- which(AICs==min(AICs))
  bic_n <- which(BICs==min(BICs))
  #print (aic_n)
  #print (bic_n)
  
  return(models[[aic_n]])
}


plot_clustering <- function(clusters, fname) {
  
  
  plot_path <- fname
  pdf(plot_path, width = 15, height = 15, useDingbats = F)
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
      #coefs <- paste(all_genes_lm[row.names(clusters)[cnt],], collapse = ",")
      merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
      best_lm <- plot_scatter_and_gam2(merged_tab, row.names(clusters)[cnt], T, nrow(res_dif_all), T, label = paste0('cluster ', cur_clst))
    }
    cnt = cnt + 1
  }
  dev.off()
}

do_clustering <- function(dif_gr, num_k, method, distance = ""){
  if (method == "kmeans"){
    # kmeans
    km <- kmeans(dif_gr, num_k, nstart = 50)
    cl <- km$cluster
  } else if (method == "pam"){
    # Pam
    if (distance == "cor"){
      cor_dif_gr <- cor(t(dif_gr))
      pam_res <- cluster::pam(cor2dist(cor_dif_gr), num_k, diss = T)
      cl <- pam_res$clustering
    } else if (distance == "eucl"){
      pam_res <- cluster::pam(scale(dif_gr), num_k, diss = F, metric = "euclidean")
      cl <- pam_res$clustering
    }
    else {
      pam_res <- cluster::pam(dif_gr, num_k, diss = F, metric = "euclidean")
      cl <- pam_res$clustering
    }
    #clusters <- as.data.frame(sort(pam_res$clustering))
    #plot_clustering(clusters, paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_gr3-_", di, "_", num_k, ".diff.sorted.cor.pdf"))
    return (cl)
  }
  
}

find_optimal_k <- function(dif_gr, method){
  if (method == "pam"){
    n_clust <- fviz_nbclust(dif_gr, pam, method = "silhouette", k.max = min(10, nrow(dif_gr) - 1))
  } else if (method == "kmeans"){
    n_clust <- fviz_nbclust(dif_gr, kmeans, method = "silhouette", k.max = min(10, nrow(dif_gr) - 1))
  } else {
    print("Wrong method! Exiting")
    return(NULL)
  }
  #n_clust
  n_clust<-n_clust$data
  max_cluster<-as.numeric(n_clust$clusters[which.max(n_clust$y)])
  if (max_cluster < 5){
    num_k <- max_cluster
  } else {
    num_k = 2
  }
  return (num_k)
}

reduce_pnts <- function(bins){
  #delta = (75-20)/300
  if (all(bins[,3] == c(0, 1, 1, 0))) bins[,3] <- c(0, 0, 1, 0)
  if (all(bins[,3] == c(1, 1, 1, 0))) bins[,3] <- c(0, 0, 1, 0)
  if (all(bins[,3] == c(0, 0, 1, 1))) bins[,3] <- c(0, 0, 1, 0)
  if (all(bins[,3] == c(1, 1, 0, 0))) bins[,3] <- c(0, 1, 0, 0)
  if (all(bins[,3] == c(0, 1, 1, 1))) bins[,3] <- c(0, 1, 0, 0)
  return(bins)
}

find_intersection_points <- function(res_dif){
  bins <- matrix(ncol = 3, nrow = 4)
  colnames(bins) <- c("start", "end", "present")
  bins[1,1] <- 0
  bins[1,2] <- 35
  bins[2,1] <- 35
  bins[2,2] <- 45
  bins[3,1] <- 45
  bins[3,2] <- 60
  bins[4,1] <- 60
  bins[4,2] <- 100
  bins[,3] <- FALSE
  #bins[,3] <- paste0(bins[,1], "_", bins[,2])
  #bins <- as.numeric(bins)
  age <- seq(20, 75, length = 300)
  #delta <- (max(res_dif) - min(res_dif))/300
  #delta <-  max((res_dif - lag(res_dif)), na.rm = T)
  delta = 0.005
  pnts <- as.numeric(age[which(abs(res_dif) < delta)])
  print(paste(row.names(res_dif), pnts))
  res_bins <- c()
  if (any(pnts > bins[1,1] & pnts < bins[1,2])) bins[1,3] = T
  if (any(pnts > bins[2,1] & pnts < bins[2,2])) bins[2,3] = T
  if (any(pnts > bins[3,1] & pnts < bins[3,2])) bins[3,3] = T
  if (any(pnts > bins[4,1] & pnts < bins[4,2])) bins[4,3] = T
  # TODO: take the min of the consecutive!!
  bins <- reduce_pnts(bins)
  
  plot_intersection_rect(res_dif, bins, delta)
  
  return(which(bins[,3] == T))
  #return (bins)
}

plot_intersection_rect <- function(res_dif, bins, delta){
  age <- seq(20, 75, length = 300)
  plot(age, res_dif, type = 'l', main =  row.names(res_dif))
  abline(h=0, lty = 1)
  abline(h=-1*delta, lty = 2)
  abline(h=delta, lty = 2)
  for (b in 1:nrow(bins)){
    if (bins[b,3]){
      print (b)
      rect(bins[b,1],-5,bins[b,2],5, col= rgb(0,0,1,alpha=0.2))
    }
  }
}
############################################################################

#x <- res_dif_all[,c("age","TFPI")]
#x <- res_dif_all[,c("age","FKBP1B")]

res_dif_all_genes <- read.delim("expression_selected_res_table.diff.txt", header = T, sep = "\t", as.is = T, check.names = F)
rand = sample(3:ncol(res_dif_all_genes), 10, replace=F)
res_dif_all <- res_dif_all_genes[, c(2,rand)]

#res_dif_all <- read.delim("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.txt", header = T, sep = "\t", as.is = T, check.names = F)
#res_dif_all <- subset(res_dif_all, select = -1)
t_res_dif_all <- as.data.frame(t(res_dif_all))

grouped_res_dif_all <- data.frame(matrix(nrow = ncol(res_dif_all) - 1, ncol = 4))
colnames(grouped_res_dif_all) <- c("degree", "cluster", "indices", "intersection")
row.names(grouped_res_dif_all) <- colnames(res_dif_all)[-1]

for (g in colnames(res_dif_all)[-1]){
  #print(g)
  d <- get_poly_degree(g)
  grouped_res_dif_all[g, "degree"] <- d
  
  int_pnts <- find_intersection_points(t_res_dif_all[g,])
  print(int_pnts)
  grouped_res_dif_all[g, "intersection"] <- paste(int_pnts, collapse = ",")
}



grouped_dif <- grouped_res_dif_all[order(grouped_res_dif_all$degree),]

grouped_dif$indices <- 0
gene_conv_table <- gene_table[match(row.names(grouped_dif), gene_table[,2]),]
all(gene_conv_table[,2] == row.names(grouped_dif))
grouped_dif$indices <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)

grouped_dif[grouped_dif$degree == 2, "degree"] = 1 
grouped_dif[grouped_dif$degree > 3, "degree"] = 4 # set all degrees > 3 to 4, change later


cl_method = "pam"
cl_dist = "cor"
for (d in unique(grouped_dif$degree)){
  dif_gr <- t_res_dif_all[row.names(grouped_dif[grouped_dif$degree == d,]),]
  
  num_k <- find_optimal_k(dif_gr, method = cl_method)
  clust <- do_clustering(dif_gr, num_k, method = cl_method, distance = cl_dist)
  grouped_dif[grouped_dif$degree == d, "cluster"] = paste0(d, "_", clust)
}

grouped_dif$cluster <- paste0(grouped_dif$cluster, "_", grouped_dif$intersection)
grouped_dif <- grouped_dif[order(grouped_dif$cluster),]
plot_clustering(grouped_dif, paste0("../anova+gam/rnaseq_plots/ageDEgenes/new/LLD_expression_random_genes_clustering_intersection_", cl_method, "_", cl_dist, "_2.pdf"))



