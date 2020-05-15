library(SimilarityMeasures)
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")

setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/all_LL")

genes_diff_all <- as.data.frame(t(read.delim("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
pheno_diff_all <- as.data.frame(t(read.delim("Laboratory_assessment_Blood_1A.fasting.datdiff.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

pheno_name = "AF"
for (g in row.names(genes_diff_all)){
  c <- cor(as.numeric(genes_diff_all[g,]), as.numeric(pheno_diff_all[pheno_name,]))
  if(c  > 0.90){
    print(paste(g, c))
  }
}



############################################
########## Compare per gender trends #######
############################################


traits0 <- as.data.frame(t(read.table(gzfile("../data/LLD_genelevel_readcounts.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.sign_inter.txt.gz"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
gte <- read.table("../data/gte_all.txt", sep = "\t", as.is = T, check.names = F)
gte_m <- gte[match(row.names(traits0), gte[,2], nomatch = 0 ),]
traits <- traits0[match(gte_m[,2], row.names(traits0), nomatch = 0 ),]
all(row.names(traits) == gte_m[,2])
row.names(traits) = gte_m[,1]

pheno0 <- as.data.frame(t(read.table("../data/age_gender_smk_contrac_cell_counts.cc_outliers_na.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
pheno <- na.omit(pheno0)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

gene_table_path <- "../data/geneid_to_gene_proteincoding_mainchr.txt"
gene_table <- read.table(gene_table_path, header = T, sep = "\t", as.is = T, check.names = F)

##############################
#############
####

genes <- c('CD248', 'ELANE', 'PTCRA', 'SELP', 'GPA33', 'TMEM40')
phenos <- c('AF', 'CHO', 'FOS')

pdf("test_distances2.pdf", width = 15, height = 15)
par(mfrow=c(3,2)) 
for (pheno_name in phenos){
  print(pheno_name)
  #pheno_name = 'AF'
  pdat_pheno <- read.delim(paste0("tmp_fitted_", pheno_name, ".txt"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
  pdat_pheno$pred_scaled <- NA
  pdat_pheno <- pdat_pheno[pdat_pheno$age > 30 & pdat_pheno$age < 70,]
  #pdat_pheno[pdat_pheno$gender_F1M2 == 1,"pred_scaled"] <- scale(pdat_pheno[pdat_pheno$gender_F1M2 == 1,'pred'])
  #pdat_pheno[pdat_pheno$gender_F1M2 == 2,"pred_scaled"] <- scale(pdat_pheno[pdat_pheno$gender_F1M2 == 2,'pred'])
  pdat_pheno$pred_scaled <- scale(pdat_pheno$pred)
  
  for (g in genes){
    print(g)
    g_id <- gene_table[gene_table$`Gene name` == g,1]
    idx = match(g_id,colnames(traits_m))
    
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = F)
    merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, g, correctForCellCounts = T, n_points = 300, make_plots = F, gam_family = gaussian(), label = '')
    pdat_g <- res_dif_lst[["pdat"]]
    pdat_g <- pdat_g[pdat_g$age > 30 & pdat_g$age < 70,]
    pdat_g$pred_scaled <- scale(pdat_g$pred)
    
    w1 <- pdat_pheno[pdat_pheno$gender_F1M2 == 1, "pred_scaled"]
    w2 <- pdat_g[pdat_g$gender_F1M2 == 1, "pred_scaled"]
    
    m1 <- pdat_pheno[pdat_pheno$gender_F1M2 == 2, "pred_scaled"]
    m2 <- pdat_g[pdat_g$gender_F1M2 == 2, "pred_scaled"]
    
    cor_w <- cor(w1,w2)
    cor_m <- cor(m1,m2)
    
    #fre_w <- Frechet(w1,w2)
    #fre_m <- Frechet(m1,m2)
    
    euc_w <- eucl_dist(w1,w2)
    euc_m <- eucl_dist(m1,m2)
    
    pl_title <- paste0("cor_w = ", format(cor_w, digits = 2), ", cor_m = ", format(cor_m, digits = 2),
                       "\nfre_w = ", format(fre_w, digits = 4), ", fre_m = ", format(fre_m, digits = 4),
                       "\neuc_w = ", format(euc_w, digits = 4), ", euc_m = ", format(euc_m, digits = 4))
    
    plot(1, type="n", main = pl_title, xlab="Age", ylab=paste0(g, " - ", pheno_name), xlim=c(20, 75), ylim=c(min(pdat_g$pred, pdat_pheno$pred_scaled), max(pdat_g$pred, pdat_pheno$pred_scaled)))
    
    levs <- levels(merged_tab$gender_F1M2)
    cols = c("indianred1", "dodgerblue1")
    for (l in seq_along(levs)) {
      dd1 <- pdat_g[pdat_g$gender_F1M2 == levs[l],]
      lines(pred_scaled ~ age, data = dd1, col = cols[l], lwd = 2)
      
      dd2 <- pdat_pheno[pdat_pheno$gender_F1M2 == levs[l],]
      lines(pred_scaled ~ age, data = dd2, col = cols[[l]], lwd = 2, lty = 2)
      
    }
  }
}
dev.off()

eucl_dist <- function(v1, v2){
  sqrt(sum((v1-v2)^2))
}


# fre < 1.2
# euc < 10
# cor > 0.8 ?


#############################
#######################
###################
n_points=300
genes_fitted_all0 <- read.delim("../anova+gam/rnaseq_plots/all_signif_interactions.fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno_fitted_all0 <- read.delim("Laboratory_assessment_Blood_1A.fasting.dat.fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

pheno_fitted_all <- as.data.frame(scale(pheno_fitted_all0))
genes_fitted_all <- as.data.frame(scale(genes_fitted_all0))

genes_fitted_all0$gender_F1M2 <- genes_fitted_all$gender_F1M2 <- c(rep(1,300), rep(2,300))
pheno_fitted_all0$gender_F1M2 <- pheno_fitted_all$gender_F1M2 <- c(rep(1,300), rep(2,300))

pheno_fitted_all0$age <- pheno_fitted_all$age <- c(seq(20, 75, length = n_points), seq(20, 75, length = n_points))
genes_fitted_all0$age <- genes_fitted_all$age <- c(seq(20, 75, length = n_points), seq(20, 75, length = n_points))

pheno_fitted_all <- pheno_fitted_all[,c(ncol(pheno_fitted_all),1:(ncol(pheno_fitted_all)-1))]
genes_fitted_all <- genes_fitted_all[,c(ncol(genes_fitted_all),1:(ncol(genes_fitted_all)-1))]

fre_w <- fre_m <- NA

pheno_fitted_all <- pheno_fitted_all[pheno_fitted_all$age >30 & pheno_fitted_all$age < 65,]
genes_fitted_all <- genes_fitted_all[genes_fitted_all$age >30 & genes_fitted_all$age < 65,]
age0 <- seq(20, 75, length = n_points)
age <- age0[age0 > 30 & age0 < 65]

min_cor_w = 1
min_cor_m =1
min_euc_w  <- min_euc_m <- 1 

pdf("best_genes_per_pheno.pdf", width = 15, height = 15)

for (pheno_name in colnames(pheno_fitted_all)[3:ncol(pheno_fitted_all)]){
  print(pheno_name)
  par(mfrow=c(5,4))
  w0 <- pheno_fitted_all0[pheno_fitted_all0$gender_F1M2 == 1, pheno_name]
  m0 <- pheno_fitted_all0[pheno_fitted_all0$gender_F1M2 == 2, pheno_name]
  
  plot(1, type="n", main = pheno_name, xlab="Age", ylab=pheno_name, xlim=c(20, 75), ylim=c(min(w0,m0) - 3*sd(w0), max(w0,m0) + 3*sd(w0)))
  
  levs <- c(1,2)
  cols = c("indianred1", "dodgerblue1")
  lines(w0 ~ age0, col = cols[1], lwd = 2)
  lines(m0 ~ age0, col = cols[2], lwd = 2)
  
  

  #pheno_name = "CHO"
  for (g in colnames(genes_fitted_all)[3:ncol(genes_fitted_all)]){
    w1 <- pheno_fitted_all[pheno_fitted_all$gender_F1M2 == 1, pheno_name]
    w2 <- genes_fitted_all[genes_fitted_all$gender_F1M2 == 1, g]
    
    m1 <- pheno_fitted_all[pheno_fitted_all$gender_F1M2 == 2, pheno_name]
    m2 <- genes_fitted_all[genes_fitted_all$gender_F1M2 == 2, g]
    
    cor_w <- cor(w1,w2)
    cor_m <- cor(m1,m2)
  
    euc_w <- eucl_dist(w1,w2)
    euc_m <- eucl_dist(m1,m2)
    
    if (!is.na(cor_m) & !is.na(cor_w)) {
      #if (euc_w < min_euc_w) min_euc_w <- euc_w
      #if (euc_m < min_euc_m) min_euc_m <- euc_m
      #if (cor_w < min_cor_w) min_cor_w <- cor_w
      #if (cor_m < min_cor_m) min_cor_m <- cor_m
      
      
        #fre_w <- Frechet(as.matrix(w1),as.matrix(w2))
        #fre_m <- Frechet(as.matrix(m1),as.matrix(m2))
      
      
      if (cor_w > 0.6 & cor_m > 0.6 & euc_m < 6 & euc_w < 6){
        # plot the original expression plot
        g_id <- gene_table[gene_table$`Gene name` == g,1]
        idx = match(g_id,colnames(traits_m))
        if (length(idx > 1)) idx <- idx[!is.na(idx)]
        #print(g)
        
        merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = F)
        merged_tab$gender_F1M2 <- as.factor(merged_tab$gender_F1M2)
        res_dif_lst <- plot_scatter_and_gam2(merged_tab, g, correctForCellCounts = T, n_points = 300, make_plots = T, gam_family = gaussian(), label = '')
        
        # plot pheno and expression together
        pl_title <- paste0("cor_w = ", format(cor_w, digits = 2), ", cor_m = ", format(cor_m, digits = 2),
                           "\neuc_w = ", format(euc_w, digits = 4), ", euc_m = ", format(euc_m, digits = 4))
        
        plot(1, type="n", main = pl_title, xlab="Age", ylab=paste0(g, " - ", pheno_name), xlim=c(20, 75), ylim=c(min(w1,w2,m1,m2), max(w1,w2,m1,m2)))
        
        levs <- c(1,2)
        cols = c("indianred1", "dodgerblue1")
        lines(w1 ~ age, col = cols[1], lwd = 2)
        lines(m1 ~ age, col = cols[2], lwd = 2)
        lines(w2 ~ age, col = cols[1], lwd = 2, lty = 2)
        lines(m2 ~ age, col = cols[2], lwd = 2, lty = 2)
        #print(g)
      
      }
    }
  }
}
dev.off()



