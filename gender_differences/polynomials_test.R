
plot_scatter_gam_lm <- function(merged_tab, trait_name, n_points, make_plots, label = ""){
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  women <- merged_tab$gender_F1M2 == 1
  lm1 <- findBestModel(merged_tab)
  
  coefs = fill_coefficients_bin(lm1, coef_names)
  lm1_formula <- paste0(formatC(lm1$coefficients, format = 'f', digits = 2), "*", names(lm1$coefficients) , collapse = " + ")
  
  if (make_plots){
    m1 <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
    m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
               data = merged_tab, method = 'REML')
    
    pdat_gam <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                             ba = mean(ba), eo = mean(eo), er = mean(er), gr = mean(gr), 
                                             ly = mean(ly),  mo = mean(mo), tr = mean(tr), gender_F1M2 = c('1', '2')))
    pdat_gam <- transform(pdat_gam, pred = predict(m1, newdata = pdat_gam, type = "response"))
    
    
    
    
    pdat_lm <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                            gender_F1M2 = c('1', '2'), ba = mean(ba), eo = mean(eo), er = mean(er), gr = mean(gr), 
                                            ly = mean(ly),  mo = mean(mo), tr = mean(tr)))
    pdat_lm <- transform(pdat_lm, pred = predict(lm1, newdata = pdat_lm, type = "response"))
    
    ## draw base plot
    ylims <- with(merged_tab, range(phenotype))
    palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
    plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
         main = paste0(pheno_name, "\n", paste(fill_coefficients_bin(lm1, coef_names), collapse = ","), "\n", label), cex = 0.6, cex.main = 0.5, xlab = "age", ylab = pheno_name)
    
    levs <- levels(merged_tab$gender_F1M2)
    cols = c("indianred1", "dodgerblue1")
    
    ## add the fitted lines
    for (l in seq_along(levs)) {
      #dd <- pdat_gam[pdat_gam$gender_F1M2 == levs[l],]
      #lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
      
      dd2 <- pdat_lm[pdat_lm$gender_F1M2 == levs[l],]
      lines(pred ~ age, data = dd2, col = cols[[l]], lwd = 2)
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
    if ("gender_F1M22:age" %in% names(lm_fit$coefficients) |  "gender_F1M22:I(age^2)" %in% names(lm_fit$coefficients) | 
                 "gender_F1M22:I(age^3)" %in% names(lm_fit$coefficients)){
      res_coef[cnt, "InteractionPresent"] == T
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
      res_coef[cnt] <- sign(lm_fit$coefficients[c])
    }
    if ("gender_F1M22:age" %in% names(lm_fit$coefficients) | "gender_F1M22:I(age^2)" %in% names(lm_fit$coefficients) | 
        "gender_F1M22:I(age^3)" %in% names(lm_fit$coefficients)){
      res_coef[cnt, "InteractionPresent"] == T
    }
    cnt = cnt + 1
  }
  return(res_coef)
}


findBestModel <- function(merged_tab){
  lm12 <- lm(phenotype ~ gender_F1M2 + age + I(age^3) + gender_F1M2:(age + I(age^3)) + ba + eo + er + gr + 
               ly + mo + tr, data = merged_tab) 
  lm11 <- lm(phenotype ~ gender_F1M2 + age + I(age^3) + ba + eo + er + gr + 
               ly + mo + tr, data = merged_tab) 
  
  lm10 <- lm(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)) + ba + eo + er + gr + 
                                                                                  ly + mo + tr, data = merged_tab) 
  lm9 <- lm(phenotype ~ gender_F1M2 + age + I(age^2) + gender_F1M2:(age + I(age^2)) + ba + eo + er + gr + 
                                                                      ly + mo + tr, data = merged_tab)
  lm8 <- lm(phenotype ~ gender_F1M2 + age + gender_F1M2:age + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  
  lm7 <- lm(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  lm6 <- lm(phenotype ~ gender_F1M2 + age + I(age^2) + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  lm5 <- lm(phenotype ~ gender_F1M2 + age + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  
  lm4 <- lm(phenotype ~ age + I(age^2) + I(age^3) + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  lm3 <- lm(phenotype ~ age + I(age^2) + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  lm2 <- lm(phenotype ~ age + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  
  lm1 <- lm(phenotype ~ gender_F1M2 + ba + eo + er + gr + 
              ly + mo + tr, data = merged_tab)
  
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


coef_names <- c( "gender_F1M22", "age", "I(age^2)", "I(age^3)", "gender_F1M22:age", "gender_F1M22:I(age^2)", 
                 "gender_F1M22:I(age^3)", "InteractionPresent", "IntersectionPoints")
all_genes_lm <- as.data.frame(matrix(ncol = length(coef_names), nrow = nrow(t_res_dif_all)))
row.names(all_genes_lm) <- row.names(t_res_dif_all)
colnames(all_genes_lm) <- coef_names

cnt = 1
#pdf("../anova+gam/rnaseq_plots/ageDEgenes/gam_vs_lm_aic.pdf", width = 15, height = 15)
#par(mfrow=c(5,5))
for (g in row.names(t_res_dif_all)){
  #g <- "GPA33"
  print(g)
  g_id <- gene_table[gene_table[,2] == g,1]
  idx <- which(colnames(traits_m) == g_id)
  n_points = 300
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
  
  best_lm <- plot_scatter_gam_lm(merged_tab, g, 300, F)
  all_genes_lm[g,] <- fill_coefficients_bin(best_lm, coef_names)
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
num_k = 5
di = "kmeans"
km <- kmeans(all_genes_lm, num_k, nstart = 25)
clusters <- as.data.frame(sort(km$cluster))
plot_clustering(clusters, paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_corrected_bin_", di, "_", num_k, ".diff.sorted.cor.pdf"))



# Pam
num_k = 6
di = "pam"

pam_res <- cluster::pam(all_genes_lm, num_k)
clusters <- as.data.frame(sort(pam_res$clustering))
plot_clustering(clusters, paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_corrected_bin_", di, "_", num_k, ".diff.sorted.cor.pdf"), all_genes_lm)

# Weighted k-means
di = "weighted_kmeans"
num_k = 7
cl <- cclust(all_genes_lm0, k=num_k, save.data=TRUE,weights =c(1,1, 0.1, 0.1, 1, 0.1, 0.1),method="hardcl")
clusters <- as.data.frame(sort(cl@cluster))
plot_clustering(clusters, paste0("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result_", di, "_", num_k, ".diff.sorted.cor.pdf"))


plot_clustering <- function(clusters, fname, all_genes_lm) {
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
      #coefs <- paste(all_genes_lm[row.names(clusters)[cnt],], collapse = ",")
      merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
      best_lm <- plot_scatter_gam_lm(merged_tab, row.names(clusters)[cnt], nrow(res_dif_all), T, label = paste0('cluster ', cur_clst))
    }
    cnt = cnt + 1
  }
  dev.off()
}




plot.new()
plot(1, type="n", xlab=NULL, ylab="",  xlim = c(1,7), ylim=c(min(all_genes_lm), max(all_genes_lm)), cex.axis = 0.01)
axis(1, at=1:7, labels=colnames(all_genes_lm), cex.axis=0.6, las = 2)

clusters <- as.data.frame(km$cluster)

cols = brewer.pal(n = num_k, name = "Set1")
for (i in 1:nrow(all_genes_lm)){
  g = row.names(all_genes_lm)[i]
  lines(x = 1:ncol(all_genes_lm), y = as.numeric(all_genes_lm[i,]), type = 'l', col = cols[clusters[g,1]])
}
