
plot_scatter_gam_lm <- function(merged_tab, trait_name, n_points, make_plots, label = ""){
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  women <- merged_tab$gender_F1M2 == 1
  lm1 <- findBestModel(merged_tab)
   
  if (make_plots){
    m1 <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
    m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
               data = merged_tab, method = 'REML')
    
    pdat_gam <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                             gender_F1M2 = c('1', '2')))
    pdat_gam <- transform(pdat_gam, pred = predict(m1, newdata = pdat_gam, type = "response"))
    
    
    
    
    pdat_lm <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                            gender_F1M2 = c('1', '2')))
    pdat_lm <- transform(pdat_lm, pred = predict(lm1, newdata = pdat_lm, type = "response"))
    
    ## draw base plot
    ylims <- with(merged_tab, range(phenotype))
    palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
    plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
         main = paste0(pheno_name, "\n", lm1$call[2], "\n", label), cex = 0.6, cex.main = 0.5, xlab = "age", ylab = pheno_name)
    
    levs <- levels(merged_tab$gender_F1M2)
    cols = c("indianred1", "dodgerblue1")
    
    ## add the fitted lines
    for (l in seq_along(levs)) {
      dd <- pdat_gam[pdat_gam$gender_F1M2 == levs[l],]
      lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
      
      dd2 <- pdat_lm[pdat_lm$gender_F1M2 == levs[l],]
      lines(pred ~ age, data = dd2, col = cols[[l]], lwd = 2, lty = 2)
    }
  }
  return(lm1)
  
}

fill_coefficients <- function(lm_fit, coef_names){
  res_coef <- rep(0, length(coef_names))
  cnt = 1
  for (c in coef_names){
    #print(c)
    if (c %in% names(lm_fit$coefficients)){
      res_coef[cnt] <- lm_fit$coefficients[c]
      
    }
    cnt = cnt + 1
  }
  return(res_coef)
}

fill_coefficients_bin <- function(lm_fit, coef_names){
  res_coef <- rep(0, length(coef_names))
  cnt = 1
  for (c in coef_names){
    #print(c)
    if (c %in% names(lm_fit$coefficients)){
      res_coef[cnt] <- lm_fit$coefficients[c]
    }
    cnt = cnt + 1
  }
  return(res_coef)
}


findBestModel <- function(merged_tab){
  lm12 <- lm(phenotype ~ gender_F1M2 + age + I(age^3) + gender_F1M2:(age + I(age^3)), data = merged_tab) 
  lm11 <- lm(phenotype ~ gender_F1M2 + age + I(age^3), data = merged_tab) 
  
  lm10 <- lm(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)), data = merged_tab) 
  lm9 <- lm(phenotype ~ gender_F1M2 + age + I(age^2) + gender_F1M2:(age + I(age^2)), data = merged_tab)
  lm8 <- lm(phenotype ~ gender_F1M2 + age + gender_F1M2:age, data = merged_tab)
  
  lm7 <- lm(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3), data = merged_tab)
  lm6 <- lm(phenotype ~ gender_F1M2 + age + I(age^2), data = merged_tab)
  lm5 <- lm(phenotype ~ gender_F1M2 + age, data = merged_tab)
  
  lm4 <- lm(phenotype ~ age + I(age^2) + I(age^3), data = merged_tab)
  lm3 <- lm(phenotype ~ age + I(age^2), data = merged_tab)
  lm2 <- lm(phenotype ~ age, data = merged_tab)
  
  lm1 <- lm(phenotype ~ gender_F1M2, data = merged_tab)
  
  models <- list(lm1, lm2, lm3, lm4, lm5, lm6, lm7, lm8, lm9, lm10, lm11, lm12)
  AICs = NULL
  BICs = NULL
  for(i in 1:length(models)){
    AICs[i] = AIC(models[[i]])
    BICs[i] = BIC(models[[i]])
  }
  
  #see which models were chosen as best by both methods
  aic_n <- which(AICs==min(AICs))
  bic_n <- which(BICs==min(BICs))
  
  return(models[[aic_n]])
}


coef_names <- c( "gender_F1M22", "age", "I(age^2)", "I(age^3)", "gender_F1M22:age", "gender_F1M22:I(age^2)", "gender_F1M22:I(age^3)")
all_genes_lm <- as.data.frame(matrix(ncol = length(coef_names), nrow = nrow(t_res_dif_all)))
row.names(all_genes_lm) <- row.names(t_res_dif_all)
colnames(all_genes_lm) <- coef_names

cnt = 1
#pdf("../anova+gam/rnaseq_plots/ageDEgenes/gam_vs_lm_aic.pdf", width = 15, height = 15)
par(mfrow=c(5,5))
for (g in row.names(t_res_dif_all)){
  #g <- "GPA33"
  print(g)
  g_id <- gene_table[gene_table[,2] == g,1]
  idx <- which(colnames(traits_m) == g_id)
  n_points = 300
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  
  best_lm <- plot_scatter_gam_lm(merged_tab, g, 300, F)
  all_genes_lm[g,] <- fill_coefficients(best_lm, coef_names)
  cnt = cnt + 1
}
#dev.off()



####################################

# hclust
di = "euclidean"
num_k = 5
d <- dist(all_genes_lm, method = di)
C1 <- hclust(d)
plot(C1) # display dendogram
clusters <- as.data.frame(sort(cutree(C1, k = num_k)))


plot_clustering(clusters, paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_", di, "_", num_k, ".diff.sorted.cor.pdf"))


# kmeans
num_k = 4
di = "kmeans"
km <- kmeans(all_genes_lm, 5, nstart = 25)
clusters <- as.data.frame(sort(km$cluster))
plot_clustering(clusters, paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_", di, "_", num_k, ".diff.sorted.cor.pdf"))



plot_clustering <- function(clusters, fname) {
  colnames(clusters) <- c('cluster')
  clusters$indices <- 0
  gene_conv_table <- gene_table[match(row.names(clusters), gene_table[,2]),]
  all(gene_conv_table[,2] == row.names(clusters))
  clusters$indices <- match(gene_conv_table[,1], colnames(traits_m), nomatch = 0)
  
  
  plot_path <- fname
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
      res_dif_lst <- best_lm <- plot_scatter_gam_lm(merged_tab, row.names(clusters)[cnt], nrow(res_dif_all), T, label = paste0('cluster ', cur_clst))
    }
    cnt = cnt + 1
  }
  dev.off()
}