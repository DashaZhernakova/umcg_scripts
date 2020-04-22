# Microbiome
source("/groups/umcg-lld/tmp03/umcg-dzhernakova/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
traits_path <- "/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome/processed/DAG3_metaphlan_bacteria_archaea_nofiltering.txt"
st_col = 1
wd_path = "/groups/umcg-lld/tmp03/umcg-dzhernakova/gender_difs/microbiome/"
plot_basepath <- "DAG3_microbiome_gamlss_poly.pdf"

pheno_path <- "dag3_age_gender.txt"

correct_for_cellcounts = F
make_plots = T

setwd(wd_path)

# traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = "ID", sep = "\t", as.is = T, check.names = F)
traits0 <- traits0[,seq(st_col,ncol(traits0))]
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)

# Age, gender and other phenotypes
pheno0 <- read.delim(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- na.omit(pheno0)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]

pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))


traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

pheno_m$gm_mean <- gm_mean_matrix(traits_m)
# Remove bacteria with < 100 individuals
traits_m <- traits_m[ ,colSums(traits_m > 0) >= 100]
# convert to fraction
#traits_m <- traits_m/100
# make smaller than one
traits_m[traits_m == 1] <- 0.9999999

num_traits <- ncol(traits_m)
num_traits


nplotspp = 20
n_points = 300
res_dif_all <- data.frame(age = seq(20, 75, length = n_points))
res_summary <- data.frame()
res_models = list()
cnt = 1

pdf(plot_basepath, width = 15, height = 15)
par(mfrow=c(5,4)) 

indices = 1:ncol(traits_m)

cnt = 1

for (idx in indices){
  if (cnt > nplotspp & make_plots){
    cnt = 1
    par(mfrow=c(5,4))
  }
  print(idx)
  trait_id <- colnames(traits_m)[idx]
  trait_name = trait_id
  pheno_name = trait_name
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "", int_tr = F)
  
  res_dif = NULL
  
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  
  res <- test_polynomial_interaction(merged_tab, BEZI(), gm_mean_cov = T)  
  inter_pval <- res[["lr_p.val"]]
  m1 <- res[["mod"]]
  res_models[[idx]] = inter_pval
  #cnt = cnt + 1
  if (inter_pval < 0.05){
    pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                         gender_F1M2 = c('1', '2'), gm_mean = mean(gm_mean)))
    pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
    pdat_dif <- expand.grid(age = seq(20, 75, length = n_points),
                            gender_F1M2 = c('1', '2'))
    res_dif <- simple_diff(pdat)
    
    ylims = c(min(0, merged_tab[,1], pdat$pred), max(merged_tab[,1], pdat$pred))
    
    
    ylabel <- pheno_name
    if (nchar(pheno_name) > 40){
      spl <- unlist(strsplit(pheno_name, "\\|"))
      ylabel <- spl[length(spl)]
      cex_main <- 0.8
      if (nchar(pheno_name) > 50) cex_main <- 0.7
      if (nchar(pheno_name) > 60) cex_main <- 0.6
      
    }
    
    
    ## draw base plot
    palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
    plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
         main = paste0(pheno_name, ' ', label, "\nLR.test interaction p = ", format(inter_pval, digits = 3)), 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, ylim = ylims)
    
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
}

dev.off()


gm_mean_matrix = function(traits_m){
  interest_matrix <- traits_m[,grepl("s__", colnames(traits_m))]
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(interest_matrix, 1, gm_mean)
}
