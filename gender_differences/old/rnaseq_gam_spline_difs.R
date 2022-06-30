args <- commandArgs(trailingOnly = TRUE)
library(RColorBrewer)
library('dplyr')
library('mgcv')

rm_na_outliers <- function(traits_m, pheno_m, idx){
  traits_na.rm <- traits_m[!is.na(traits_m[,idx]),idx]
  pheno_na.rm <- pheno_m[!is.na(traits_m[,idx]),]
  women <- pheno_na.rm$gender_F1M2 == 1
  
  #merged_tab <- cbind(traits_na.rm, pheno_na.rm[,c("age", "gender_F1M2", "contrac")])
  #colnames(merged_tab) <- c(colnames(traits_m)[idx], "age", "gender_F1M2", "contrac")
  all(row.names(traits_na.rm) == row.names(pheno_na.rm))
  
  merged_tab <- cbind(traits_na.rm, pheno_na.rm)
  colnames(merged_tab) <- c(colnames(traits_m)[idx], colnames(pheno_na.rm))
  
  row.names(merged_tab) <- row.names(pheno_na.rm)
  
  # remove outliers ( skip this if outcome is a factor )
  if (length(unique(merged_tab[,1])) > 3){
    w <- merged_tab[women,]
    m <- merged_tab[!women,]
    
    # Zscore < 3
    #tab_nooutliers <- rbind(w[abs(w[,1] - mean(w[,1]))/sd(w[,1]) < 3,], m[abs(m[,1] - mean(m[,1]))/sd(m[,1]) < 3,])
    
    # less than 1.5 IQR from the 1st and 3rd quantiles
    mq1 <- quantile(m[,1], probs = 0.25)
    mq3 <- quantile(m[,1], probs = 0.75)
    miqr <- mq3 - mq1
    m_clean <- m[m[,1] < mq3 + 1.5*miqr & m[,1] > mq1 - 1.5*miqr,]
    
    wq1 <- quantile(w[,1], probs = 0.25)
    wq3 <- quantile(w[,1], probs = 0.75)
    wiqr <- wq3 - wq1
    w_clean <- w[w[,1] < wq3 + 1.5*wiqr & w[,1] > wq1 - 1.5*wiqr,]
    
    tab_nooutliers <- rbind(w_clean, m_clean)
    
    
  } else {
    tab_nooutliers <- merged_tab
  }
  return(tab_nooutliers)
}

plot_scatter_and_gam2 <- function(merged_tab, pheno_name, correctForCellCounts, n_points, make_plots, min_age=20, max_age=75){
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  women <- merged_tab$gender_F1M2 == 1
  
  if (!correctForCellCounts){
    # TODO: check if binomial works
    if (length(unique(merged_tab[,1])) == 2){ # binary
      m1 <- gam(phenotype ~ gender_F1M2 + s(age, by = gender_F1M2), data = merged_tab, family=binomial(link = "logit"), method = "REML")
      m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2),
                 data = merged_tab, method = 'REML', family=binomial(link = "logit"))
    } else {
      m1 <- gam(phenotype ~ gender_F1M2 + s(age) + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
      m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2),
                 data = merged_tab, method = 'REML')
    }
    
    m_o_p <- summary(m_o)$s.pv[length(summary(m_o)$s.pv)]
    m_o_g_beta <- m_o$coefficients["ord_gender_F1M2.L"]
    m_o_g_pv <- summary(m_o)$p.pv["ord_gender_F1M2.L"]
    
    if (m_o_p > 0.05){
      return (list("dif" = NULL, "inter_p" = m_o_p,"g_beta" = m_o_g_beta, "g_pv" = m_o_g_pv))
    }
    pdat <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points),
                                         gender_F1M2 = c('1', '2')))
    pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
    
    #pdat_dif <- expand.grid(age = seq(min_age, max_age, length = n_points),
    #                    gender_F1M2 = c('1', '2'))
    #res_dif <- smooth_diff(m1, pdat_dif, '2', '1', "gender_F1M2")
    #res_dif$age = seq(min_age,max_age, length = n_points)
    res_dif <- simple_diff(pdat)
    
  } else { # Correct for cell counts
    
    m1 <- gam(phenotype ~ gender_F1M2 + ba + eo + er + gr +
                ly + mo + tr + s(age, by = gender_F1M2), data = merged_tab, method = "REML")
    
    m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + ba + eo + er + gr +
                 ly + mo + tr + s(age, by = ord_gender_F1M2),
               data = merged_tab, method = 'REML')
    m_o_p <- summary(m_o)$s.pv[length(summary(m_o)$s.pv)]
    m_o_g_beta <- m_o$coefficients["ord_gender_F1M2.L"]
    m_o_g_pv <- summary(m_o)$p.pv["ord_gender_F1M2.L"]
    
    if (m_o_p > 0.05){
      return (list("dif" = NULL, "inter_p" = m_o_p,"g_beta" = m_o_g_beta, "g_pv" = m_o_g_pv))
    }
    
    pdat <- with(merged_tab, expand.grid(age = seq(min_age, max_age, length = n_points), gender_F1M2 = c('1', '2'),
                                         ba = mean(ba), eo = mean(eo), er = mean(er), gr = mean(gr),
                                         ly = mean(ly),  mo = mean(mo), tr = mean(tr)))
    pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
    
    #pdat_dif <- expand.grid(age = seq(min_age, max_age, length = n_points), gender_F1M2 = c('1', '2'),
    #                    ba = mean(merged_tab$ba), eo = mean(merged_tab$eo), er = mean(merged_tab$er), gr = mean(merged_tab$gr),
    #                    ly = mean(merged_tab$ly),  mo = mean(merged_tab$mo), tr = mean(merged_tab$tr))
    #res_dif <- smooth_diff(m1, pdat_dif, '2', '1', "gender_F1M2")
    #res_dif$age = seq(min_age,max_age, length = n_points)
    res_dif <- simple_diff(pdat)
  }
  
  if (make_plots){
    ylims <- with(merged_tab, range(phenotype))
    ## draw base plot
    palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
    plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16,
         main = paste0(pheno_name, "\nGAM interaction p = ", format(m_o_p, digits = 3)),
         cex = 0.6, xlab = "age", ylab = pheno_name)
    
    levs <- levels(merged_tab$gender_F1M2)
    cols = c("indianred1", "dodgerblue1")
    
    ## add the fitted lines
    for (l in seq_along(levs)) {
      dd <- pdat[pdat$gender_F1M2 == levs[l],]
      lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
    }
    
    #Plot the difference
    plot(res_dif$age, res_dif$diff, type = 'l')
  }
  
  return (list("pdat" = pdat,"dif" = res_dif$diff, "inter_p" = m_o_p,"g_beta" = m_o_g_beta, "g_pv" = m_o_g_pv))
}

# Calculate a simple difference between every fitted value pair in men and women
simple_diff <- function(pdat){
  pdat <- pdat[order(pdat$gender_F1M2, pdat$age),]
  dif <- pdat[pdat$gender_F1M2 == 2, "pred"] - pdat[pdat$gender_F1M2 == 1, "pred"]
  res_dif <- data.frame(age = pdat[pdat$gender_F1M2 == 1, "age"], diff = dif)
  return(res_dif)
}

smooth_diff <- function(model, newdata, f1, f2, var, alpha = 0.05,
                        unconditional = FALSE) {
  xp <- predict(model, newdata = newdata, type = 'lpmatrix')
  c1 <- grepl(f1, colnames(xp))
  c2 <- grepl(f2, colnames(xp))
  r1 <- newdata[[var]] == f1
  r2 <- newdata[[var]] == f2
  ## difference rows of xp for data from comparison
  X <- xp[r1, ] - xp[r2, ]
  ## zero out cols of X related to splines for other lochs
  X[, ! (c1 | c2)] <- 0
  ## zero out the parametric cols
  X[, !grepl('^s\\(', colnames(xp))] <- 0
  dif <- X %*% coef(model)
  se <- sqrt(rowSums((X %*% vcov(model, unconditional = unconditional)) * X))
  crit <- qt(alpha/2, df.residual(model), lower.tail = FALSE)
  upr <- dif + (crit * se)
  lwr <- dif - (crit * se)
  res_dif <- data.frame(pair = paste(f1, f2, sep = '-'),
                        diff = dif,
                        se = se,
                        upper = upr,
                        lower = lwr)
}

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}
######################################################################################
boxy = T

if (boxy){
  wd_path <- '/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/genes_sign_interaction/'
} else {
  wd_path <- "C:/Users/Dasha/work/UMCG/data/gender_differences/omics/data/"
}
if  (length(args) > 1) {
  traits0_path <- args[1]
  out_prefix <- args[2]
} else {
  traits0_path <- "../LLD_expression/gene_read_counts_BIOS_and_LLD_passQC.only_LLD.no_zeros.TMM.CPM.tsv.gz.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"
  #traits0_path <- "gene_subset_normalized.txt"
  
  out_prefix <- "expression_signif_res_table"
}
gte_path <- "gte_all.txt"
pheno_path <- "age_gender_smk_contrac_cell_counts.cc_outliers_na.txt"
gene_table_path <- "geneid_to_gene_proteincoding.txt"
plot_basepath <- out_prefix
correct_for_cellcounts = TRUE
make_plots = F

print(paste('wd:', wd_path))
print(paste('traits:',traits0_path))
print(paste('pheno:', pheno_path))
print(paste('out_prefix:', out_prefix))
print(paste('correct for cell counts:', correct_for_cellcounts))
print(paste('make plots:', make_plots))

setwd(wd_path)

# expression
traits0 <- as.data.frame(t(read.table(gzfile(traits0_path), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
gte <- read.table(gte_path, sep = "\t", as.is = T, check.names = F)
gte_m <- gte[match(row.names(traits0), gte[,2], nomatch = 0 ),]
traits <- traits0[match(gte_m[,2], row.names(traits0), nomatch = 0 ),]
all(row.names(traits) == gte_m[,2])
row.names(traits) = gte_m[,1]
#traits <- as.data.frame(t(read.table(gzfile("/Users/dashazhernakova/Documents/UMCG/data/olink/LLDeep_expression/gene_level/LLD_genelevel_readcounts.nozeros.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

num_traits <- ncol(traits)

pheno0 <- as.data.frame(t(read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
pheno <- na.omit(pheno0)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]



gene_table <- read.table(gene_table_path, header = T, sep = "\t", as.is = T, check.names = F)


#pheno_m$gender_F1M2 <- as.factor(pheno_m$gender_F1M2)

ntraits <- ncol(traits_m)
ntraits
nplotspp <-20


if (length(args) > 3) {
  start_idx <- args[3]
  end_idx <- args[4]
  plot_path <- paste0(plot_basepath, ".", start_idx, "-", end_idx,".pdf")
  
} else {
  start_idx <- 1
  end_idx <- ncol(traits_m)
  plot_path <- paste0(plot_basepath, ".pdf")
  
}

n_points = 300
min_age = 15
max_age = 90
res_dif_all <- data.frame(age =seq(min_age,max_age, length = n_points))
res_summary <- data.frame()
res_pred_all <- data.frame(age = c(seq(min_age, max_age, length = n_points),seq(min_age,max_age, length = n_points)))
cnt = 1
if (make_plots){
  pdf(plot_path, width = 15, height = 15)
  par(mfrow=c(5,4))
}
for (idx in start_idx:end_idx){
  if (cnt > nplotspp & make_plots){
    cnt = 1
    par(mfrow=c(5,4))
  }
  print(idx)
  
  pheno_name = gene_table[gene_table[,1] == colnames(traits_m)[idx],2]
  if (length(pheno_name) > 0){ #if gene id in gene convertion table
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
    
    
    res_dif = NULL
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, pheno_name, correct_for_cellcounts, n_points, make_plots, min_age = min_age, max_age = max_age)
    
    res_summary[pheno_name,'inter_p'] = res_dif_lst[["inter_p"]]
    res_summary[pheno_name,'g_beta'] = res_dif_lst[["g_beta"]]
    res_summary[pheno_name,'g_pv'] = res_dif_lst[["g_pv"]]
    
    if (res_dif_lst[["inter_p"]] < 0.05){
      #break
      cnt <- cnt + 1
      res_dif_all[,pheno_name] <- res_dif_lst[["dif"]]
      res_pred_all[,pheno_name] <- res_dif_lst[["pdat"]]$pred
    }
  }
  
}
if (make_plots){
  dev.off()
}

write.table(res_dif_all, file = paste0(out_prefix, ".diff.txt"), sep = "\t", quote = F, col.names = NA)
write.table(res_summary, file = paste0(out_prefix, ".summary.txt"), sep = "\t", quote = F, col.names = NA)
write.table(res_pred_all, file = paste0(out_prefix, ".fitted.txt"), sep = "\t", quote = F, col.names = NA)