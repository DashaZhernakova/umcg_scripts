library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(mgcv)
library(cluster)
library(factoextra)
library(psych)

source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")

# Get the degree of the best fitted polynomial (for the difference trend)
get_poly_degree <- function (res_dif_all, g){
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


# Choose the best model for determining the polynomial for the difference trend. Needs to be redone!!!
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


# Plot traits and differences in the order determined by clustering
plot_clustering <- function(clusters, fname, scale_traits = F) {
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
      if (scale_traits){
        merged_tab <- rm_na_outliers(scale(traits_m), pheno_m, idx)
      } else {
        merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
      }
      
      best_lm <- plot_scatter_and_gam2(merged_tab, row.names(clusters)[cnt], T, nrow(res_dif_all), T, label = paste0('cluster ', cur_clst))
    }
    cnt = cnt + 1
  }
  dev.off()
}

# Plot traits and differences in the order determined by clustering usng ggplot (plots are saved in the first run together with difference calculation)
# Doesn't work!
plot_clustering_ggplot <- function(clusters, fname, plot_list) {
  plot_path <- fname
  pdf(plot_path, width = 15, height = 15, useDingbats = F)
  #par(mfrow=c(5,4))
  
  cur_clst <- 1
  cnt <- 1
  cur_page <- list()
  
  for (idx in clusters$indices){
    print (cnt)
    if (clusters[cnt,'cluster'] != cur_clst){
      cur_clst <- clusters[cnt,'cluster']
      
      #par(mfrow=c(5,4))
      gridExtra::marrangeGrob(cur_page, nrow = 5, ncol = 4)
      cur_page <- list()
      cnt = 1
    }
    if (length(cur_page) > 20){
      gridExtra::marrangeGrob(cur_page, nrow = 5, ncol = 4)
      cur_page <- list()
      cnt = 1 
    }
    if (idx != 0){
      #coefs <- paste(all_genes_lm[row.names(clusters)[cnt],], collapse = ",")
      #merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
      #best_lm <- plot_scatter_and_gam2(merged_tab, row.names(clusters)[cnt], T, nrow(res_dif_all), T, label = paste0('cluster ', cur_clst))
      
      cur_page[[cnt]] = plot_list[[idx]][[1]]
      cnt = cnt + 1
    }
    
  }
  dev.off()
}

# Perform clustering based on the difference between trends in men and in women
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

# Finds the optimal k for k-means and k-medoids clustering
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


# Find the bins where men = women (fits intersect)
find_intersection_points <- function(res_dif, degree){
  n_points = 300
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
  
  age <- seq(20, 75, length = n_points)
  
  #delta <- (max(res_dif) - min(res_dif))/300
  #delta <-  max((res_dif - lag(res_dif)), na.rm = T)
  delta = 0.005
  
  pnts <- c()
  for (i in 1:(n_points - 1)){
    l = min(res_dif[i], res_dif[i+1])
    h = max(res_dif[i], res_dif[i+1])
    #print(paste(i, age[i], l, h))
    if(h > -1*delta & l <= delta){
      pnts <- c(pnts, age[i])
    }
  }
  
  #pnts <- as.numeric(age[which(abs(res_dif) < delta)])
  #print(paste(row.names(res_dif), pnts))
  res_bins <- c()
  if (length(pnts) == 0 & degree < 3){
    pnts <- get_nearest_point(res_dif, age)
  }
  if (any(pnts > bins[1,1] & pnts < bins[1,2])) bins[1,3] = T
  if (any(pnts > bins[2,1] & pnts < bins[2,2])) bins[2,3] = T
  if (any(pnts > bins[3,1] & pnts < bins[3,2])) bins[3,3] = T
  if (any(pnts > bins[4,1] & pnts < bins[4,2])) bins[4,3] = T
  # TODO: take the min of the consecutive!!
  bins <- reduce_pnts(bins)
  
  #plot_intersection_rect(res_dif, bins, delta)
  
  return(which(bins[,3] == T))
  #return (bins)
}

# Find the point with the minimum distance between mwn and women (used when they don't intersect)
get_nearest_point <- function(res_dif, age){
  pnts <- age[which.min(abs(res_dif))]
  return(pnts)
}

# Remove the consequtive intersection bins, keep only one
reduce_pnts <- function(bins){
  #delta = (75-20)/300
  if (all(bins[,3] == c(0, 1, 1, 0))) bins[,3] <- c(0, 0, 1, 0)
  if (all(bins[,3] == c(1, 1, 1, 0))) bins[,3] <- c(0, 0, 1, 0)
  if (all(bins[,3] == c(0, 0, 1, 1))) bins[,3] <- c(0, 0, 1, 0)
  if (all(bins[,3] == c(1, 1, 0, 0))) bins[,3] <- c(0, 1, 0, 0)
  if (all(bins[,3] == c(0, 1, 1, 1))) bins[,3] <- c(0, 1, 0, 0)
  return(bins)
}

# Plot the difference and the intersection bins
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