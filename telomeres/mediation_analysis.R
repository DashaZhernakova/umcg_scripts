setwd("C:/Users/Dasha/work/UMCG/data/telomeres")
library(mediation)
library(ppcor)

traits_path <- "methylation_scaled_v3_top1k.txt.gz"
gte_path <- "../LifeLines_phenotypes/methylation_gte.txt"
gene_table_path <- "../LifeLines_phenotypes/meth_probe2gene.txt"

pheno_path1 <- "parental_age_smk_diab_final.v3.covariates.txt"
pheno_path2 <- "bioage_gender_age.txt"



traits0 <- as.data.frame(t( read.table(gzfile(traits_path), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
#row.names(traits0) <- gsub("X","", row.names(traits0))
#gte <-read.delim(gte_path,  sep = "\t", as.is = T, check.names = F, header = T )
#gte_m <- gte[match(row.names(traits0), gte[, 2], nomatch = 0), ]
#traits <- traits0[match(gte_m[, 2], row.names(traits0), nomatch = 0),]
#all(row.names(traits) == gte_m[, 2])
#row.names(traits) <- gte_m[, 1]
traits <- traits0
nrow(traits)

pheno <- read.table(pheno_path1, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0), ]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))

pheno2 <- read.table(pheno_path2, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits_m2 <- traits_m[match(row.names(pheno2), row.names(traits_m), nomatch = 0), ]
pheno_m2 <- pheno_m[match(row.names(pheno2), row.names(pheno_m), nomatch = 0), ]
pheno2_m <- pheno2[match(row.names(traits_m2), row.names(pheno2), nomatch = 0), ]
all(row.names(traits_m2) == row.names(pheno2_m))
nrow(traits_m2)

gene_table <- read.table(gene_table_path, header = T, sep = "\t", as.is = T, check.names = F)

pheno_m$gender <- as.factor(pheno_m$gender)
if ("SMK13" %in% colnames(pheno_m)) pheno_m$SMK13 <- as.factor(pheno_m$SMK13)
if ("SMK14" %in% colnames(pheno_m)) pheno_m$SMK14 <- as.factor(pheno_m$SMK14)
if ("parents_diab" %in% colnames(pheno_m)) pheno_m$parents_diab <- as.factor(pheno_m$parents_diab)
if ("father_diab" %in% colnames(pheno_m)) pheno_m$father_diab <- as.factor(pheno_m$father_diab)
if ("mother_diab" %in% colnames(pheno_m)) pheno_m$mother_diab <- as.factor(pheno_m$mother_diab)


me_parent <- read.delim("results/parental_pheno/methylation_vs_parental_bender_main_cellcount_corrected.scaled.v3.p0.05.perm_FDR.txt", header = T, sep = "\t", as.is = T, check.names = F)
me_mtl <- read.delim("results/mtl/methylation_vs_mtl_bender_main_cellcount_corrected_scaled.v3.p0.05.txt", header = T, sep = "\t", as.is = T, check.names = F)
#me_mtl <- me_mtl[me_mtl$pval < 0.05,]

res <- data.frame(matrix(ncol = 19, nrow = 0))
colnames(res) <- c("probe", "mtl", "pheno","probe-pheno", "probe-mtl", "mtl-pheno","pcor_pval.pheno-mtl", "pcor_pval.pheno-probe", "pcor_pval.probe-mtl", "pcor_correl.pheno-mtl", "pcor_correl.pheno-probe", "pcor_correl.probe-mtl", "ACME_pval", "ACME_estimate", "ADE_pval", "Prop_med", "ACME_pval_reverse", "ADE_pval_reverse", "Prop_med_reverse")

cnt <- 1
for (i in 1:nrow(me_parent)){
  probe <- me_parent[i,"pheno1"]
  pheno <- me_parent[i, "pheno2"]
  me_mtl_subs <- me_mtl[me_mtl$pheno1 == probe,]
  if (nrow(me_mtl_subs) > 0){
  for (j in 1: nrow(me_mtl_subs)){
    mtl <- me_mtl_subs[j, "pheno2"]
    
    #merged_tab <- cbind(traits_m2[,probe], pheno_m2[,c('age', 'gender', pheno)], pheno2_m[,mtl], )
    #colnames(merged_tab) <- c('probe', 'age', 'gender', 'pheno', 'mtl')
    merged_tab <- cbind(traits_m2[,probe], pheno_m2[,c('age', 'gender', "ba", "eo", "er", "gr", "ly", "mo", "tr", pheno)], pheno2_m[,mtl])
    colnames(merged_tab) <- c('probe', 'age', 'gender', "ba", "eo", "er", "gr", "ly", "mo", "tr", 'pheno', 'mtl')
    
    # print(paste(probe, mtl, pheno))
    
    res[cnt, 'probe'] <- probe
    res[cnt, 'mtl'] <- mtl
    res[cnt, 'pheno'] <- pheno
    res[cnt, "probe-pheno"] <- me_parent[i, "beta"]
    res[cnt, "probe-mtl"] <- me_mtl_subs[j, "beta"]

    #res[cnt,] <- run_mediation(merged_tab, res[cnt,])
    res[cnt,] <- run_mediation_correct_cellcounts(merged_tab, res[cnt,])
    cnt = cnt + 1
  }
  }
}
write.table(res, file = "mediation_results_main_cell_counts_v3.txt", sep = "\t", quote = F, col.names = NA)


run_mediation <- function(merged_tab, res){
  
  merged_tab <- na.omit(merged_tab)
  
  fit.totaleffect=lm(mtl ~ age + gender + pheno,merged_tab)
  res$fit.totaleffect <- summary(fit.totaleffect)$coefficients[4,4]
  res$`mtl-pheno` <- summary(fit.totaleffect)$coefficients[4,1]
  
  if (res$fit.totaleffect > 0.05)
    return (res)
  fit.mediator=lm(probe ~ age + gender + pheno,merged_tab)
  res$fit.mediator <- summary(fit.mediator)$coefficients[4,4]
  
  fit.dv=lm(mtl ~ age + gender + probe + pheno, merged_tab)
  res$fit.dv <- summary(fit.dv)$coefficients[5,4]
  
  #print(paste(p2, p1, p3))
  results = mediate(fit.mediator, fit.dv, treat='pheno', mediator='probe', boot=T)
  #print(summary(results))
  res$ACME_pval <- summary(results)$d0.p
  res$ADE_pval <- summary(results)$z0.p
  res$Prop_med <- summary(results)$n0.p
  res$ACME_estimate <- summary(results)$d0
  
  ### Opposite direction
  fit.totaleffect2=lm(probe ~ age + gender + pheno, merged_tab)
  fit.mediator2=lm(mtl ~ age + gender + pheno,merged_tab)
  fit.dv2=lm(probe ~ age + gender + mtl + pheno,merged_tab)
  results2 = mediate(fit.mediator2, fit.dv2, treat='pheno', mediator='mtl', boot=T)
  #print(summary(results2))
  res$ACME_pval_reverse <- summary(results2)$d0.p
  res$ADE_pval_reverse <- summary(results2)$z0.p
  res$Prop_med_reverse <- summary(results2)$n0.p
  
  return(res)
}

run_mediation_correct_cellcounts <- function(merged_tab, res){
  
  merged_tab <- na.omit(merged_tab)
  pcor_res <-pcor(merged_tab, method = "pearson")
  res$`pcor_pval.pheno-mtl` <- pcor_res$p.value["pheno", "mtl"]
  res$`pcor_pval.pheno-probe` <- pcor_res$p.value["pheno", "probe"]
  res$`pcor_pval.probe-mtl` <- pcor_res$p.value["mtl", "probe"]
  res$`pcor_correl.pheno-mtl` <- pcor_res$estimate["pheno", "mtl"]
  res$`pcor_correl.pheno-probe` <- pcor_res$estimate["pheno", "probe"]
  res$`pcor_correl.probe-mtl` <- pcor_res$estimate["mtl", "probe"]

  fit.totaleffect=lm(mtl ~ age + gender + pheno + ba + eo + er + gr + ly + mo + tr,merged_tab)
  res$`mtl-pheno` <- summary(fit.totaleffect)$coefficients[4,1]
  #res$fit.totaleffect <- summary(fit.totaleffect)$coefficients[4,4]
  
  if (summary(fit.totaleffect)$coefficients[4,4] > 0.05)
    return (res)
  fit.mediator=lm(probe ~ age + gender + pheno + ba + eo + er + gr + ly + mo + tr,merged_tab)
  #res$fit.mediator <- summary(fit.mediator)$coefficients[4,4]
  
  fit.dv=lm(mtl ~ age + gender + probe + pheno + ba + eo + er + gr + ly + mo + tr, merged_tab)
  #res$fit.dv <- summary(fit.dv)$coefficients[5,4]
  
  #print(paste(p2, p1, p3))
  results = mediate(fit.mediator, fit.dv, treat='pheno', mediator='probe', boot=T)
  #print(summary(results))
  res$ACME_pval <- summary(results)$d0.p
  res$ADE_pval <- summary(results)$z0.p
  res$Prop_med <- summary(results)$n0
  res$ACME_estimate <- summary(results)$d0
  
  ### Opposite direction
  fit.totaleffect2=lm(probe ~ age + gender + pheno + ba + eo + er + gr + ly + mo + tr, merged_tab)
  fit.mediator2=lm(mtl ~ age + gender + pheno + ba + eo + er + gr + ly + mo + tr, merged_tab)
  fit.dv2=lm(probe ~ age + gender + mtl + pheno + ba + eo + er + gr + ly + mo + tr, merged_tab)
  results2 = mediate(fit.mediator2, fit.dv2, treat='pheno', mediator='mtl', boot=T)
  #print(summary(results2))
  res$ACME_pval_reverse <- summary(results2)$d0.p
  res$ADE_pval_reverse <- summary(results2)$z0.p
  res$Prop_med_reverse <- summary(results2)$n0
  
  return(res)
}
for (probe in c("cg08620154", "cg18897632", "cg01566526")){
  for (mtl in c("MTL_lymp", "MTL_CD57+", "MTL_CD45+_CD20-", "MTL_gran", "MTL_CD20+", "MTL_CD45-")){

#probe = "cg18897632"
#mtl = "MTL_CD45-"
#pheno = "age_m"
merged_tab <- cbind(traits_m2[,probe], pheno_m2[,c('age', 'gender', "ba", "eo", "er", "gr", "ly", "mo", "tr", pheno)], pheno2_m[,mtl])
colnames(merged_tab) <- c('probe', 'age', 'gender', "ba", "eo", "er", "gr", "ly", "mo", "tr", 'pheno', 'mtl')


merged_tab <- na.omit(merged_tab)
pcor_res <-pcor(merged_tab, method = "pearson")
print(paste(pheno, probe, gene_table[gene_table$`HT12v4-ArrayAddress` == probe,"Gene"], mtl,
            pcor_res$p.value["pheno", "mtl"], pcor_res$p.value["pheno", "probe"], pcor_res$p.value["mtl", "probe"],
      pcor_res$estimate["pheno", "mtl"], pcor_res$estimate["pheno", "probe"], pcor_res$estimate["mtl", "probe"]))
run_mediation_correct_cellcounts(merged_tab, res)
  }
}
#
#
names= list('probe', 'pheno', 'mtl')
data_vals = c(0, "a", 0, #X2
              0, 0, 0, #X1
              "b", "c", 0) #Y
M <- matrix(nrow = 3, ncol = 3, byrow = TRUE, data = data_vals)
plotmat(M, pos = c(1, 2), curve = 0, name = names, lwd = 1,
        box.lwd = 2, cex.txt = 0.8, box.type = "square", box.prop = 1.0)
