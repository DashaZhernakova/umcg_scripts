library(compositions)

# Microbiome
traits_path <- "LLD.metaphlan.raw.tsv"
st_col = 1

plot_basepath <- "../plots_all_pheno/LLD_microbiome_gamlss.pdf"

gte_path <- "gte_all.txt"
pheno_path <- "age_gender_smk_contrac_cell_counts.cc_outliers_na.txt"
gene_table_path <- "geneid_to_gene_proteincoding_mainchr.txt"

correct_for_cellcounts = F
make_plots = T

setwd(wd_path)

# traits of interest
traits0 <- as.data.frame(t(read.table(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
traits0 <- traits0[,seq(st_col,ncol(traits0))]
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)

# Age, gender and other phenotypes
pheno0 <- as.data.frame(t(read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
pheno <- na.omit(pheno0)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]

pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))


traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]

#
# Preprocessing
#

core_matrix <- traits_m[,grepl("s__", colnames(traits_m))]

#  Log ratio transformation
traits_m <- do_clr_externalWeighting(traits_m, core_matrix)

# Remove bacteria with < 100 individuals
traits_m <- traits_m[ ,colSums(traits_m > 0) >= 100]

traits_m <- traits_m/rowSums(traits_m)


num_traits <- ncol(traits_m)
num_traits


nplotspp = 20
n_points = 300
res_dif_all <- data.frame(age = seq(20, 75, length = n_points))
res_summary <- data.frame()
res_models = list()
cnt = 1
if (make_plots){
  pdf(plot_basepath, width = 15, height = 15)
  par(mfrow=c(5,4)) 
}


indices = 1:ncol(traits_m)
indices = 150:200
cnt = 1

for (idx in indices){
  if (cnt > nplotspp & make_plots){
    cnt = 1
    par(mfrow=c(5,4))
  }
  print(idx)
  
  trait_id <- colnames(traits_m)[idx]
  trait_name = trait_id
  if (length(trait_name) > 0){ #if gene id in gene convertion table
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "", int_tr = F)
    
    #tryCatch({
      res_models[[idx]] <- plot_scatter_and_gamlss(merged_tab, trait_name, F, 300, make_plots, label = '')
    #},error=function(e) {
    #      message(paste("Fitting failed for ", idx))
   # }) 
    
    
    #res_summary[trait_id,'p1'] = res_dif_lst[["plots"]][[1]]
    #res_summary[trait_id,'p2'] = res_dif_lst[["p2"]]
    #plot_list[[idx]] = res_dif_lst[["plots"]]
    
  }
  
}

if (make_plots){
  dev.off()
  
}


do_clr_externalWeighting = function(interest_matrix, core_matrix){
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  #estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

plot_scatter_and_gamlss0 <- function(merged_tab, pheno_name, correctForCellCounts = F, n_points, make_plots, label = '', gam_family = gaussian()){
  merged_tab <- mutate(merged_tab, ord_gender_F1M2 = ordered(gender_F1M2, levels = c('1', '2')))
  merged_tab <- mutate(merged_tab, gender_F1M2 = as.factor(gender_F1M2))
  #pheno_name = gene_table[gene_table[,1] == colnames(merged_tab)[1],2]
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < 75) & (merged_tab$age >= 20),]
  women <- merged_tab$gender_F1M2 == 1
  print(is(merged_tab$gender_F1M2, "factor"))
  #m1 <- gam(phenotype ~ gender_F1M2 + s(age, bs = 'cs') + s(age, by = gender_F1M2, bs = 'cs'), data = merged_tab, method = "REML", family = gam_family)
  #m_o <- gam(phenotype ~ ord_gender_F1M2 + s(age) + s(age, by = ord_gender_F1M2), 
  #           data = merged_tab, method = 'REML', family = gam_family)
  
  m1 <- gamlss(phenotype ~ gender_F1M2 + cs(age) + pvc(age, by=gender_F1M2), data = merged_tab, family = BEZI())
  m0 <- gamlss(phenotype ~ gender_F1M2 + pb(age), data = merged_tab, family = BEZI())
  
  if (AIC(m1) < AIC(m0)){
    print("Bingo!")
    pdat <- with(merged_tab, expand.grid(age = seq(20, 75, length = n_points), 
                                         gender_F1M2 = c('1', '2')))
    pdat <- mutate(pdat, gender_F1M2 = as.factor(gender_F1M2))
    print(is(pdat$gender_F1M2, "factor"))
    pdat <- transform(pdat, pred = predict(m1, newdata = pdat, type = "response"))
    
    pdat_dif <- expand.grid(age = seq(20, 75, length = n_points),
                            gender_F1M2 = c('1', '2'))
    res_dif <- simple_diff(pdat)
    
    ylabel <- pheno_name
    if (nchar(pheno_name) > 40){
      spl <- unlist(strsplit(pheno_name, "\\|"))
      ylabel <- spl[length(spl)]
      cex_main <- 0.8
      if (nchar(pheno_name) > 50) cex_main <- 0.7
      if (nchar(pheno_name) > 60) cex_main <- 0.6
      
    }
    
    palette(c(col2transparent("indianred1", 125),col2transparent("dodgerblue1", 125)))
    plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
         main = paste0(pheno_name), 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main)
    
    levs <- levels(merged_tab$gender_F1M2)
    cols = c("indianred1", "dodgerblue1")
    
    ## add the fitted lines
    for (l in seq_along(levs)) {
      dd <- pdat[pdat$gender_F1M2 == levs[l],]
      lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
      
    }
    plot(res_dif$age, res_dif$diff, type = 'l')
    return(1)
    }
  return(0)
}