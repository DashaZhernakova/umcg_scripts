real_f <- "methylation_vs_parental_v2_bender_main_cellcount_corrected_scaled.v3.srt.p0.05.txt"
perm_f <- "methylation_vs_parental_main_cellcount_corrected_scaled.v3.permutation1-10.p0.05.srt.txt.gz" 
nperm <- 10

real_data <- read.delim(real_f, header = T, sep = "\t", as.is = T, check.names = F)
perm_data <- read.delim(gzfile(perm_f), header = T, sep = "\t", as.is = T, check.names = F)

perm_data <- perm_data[perm_data$pheno2 != 'parents_diab',]
real_data <- real_data[real_data$pheno2 != 'parents_diab',]
perm_data0 <- perm_data
real_data0 <- real_data

real_data$FDR_tmp <- NA

cur_perm_idx <- 1
for (i in 1:nrow(real_data)){
    cur_real_pval = real_data[i,'pval']
    while (perm_data[cur_perm_idx,'pval'] < cur_real_pval){
        cur_perm_idx = cur_perm_idx + 1
    }
    fdr <- ((cur_perm_idx - 1)/nperm) / (i - 1)
    #print (paste(i-1, cur_perm_idx-1, ((cur_perm_idx - 1)/nperm) / (i - 1)))
    real_data[i,'FDR_tmp'] <- fdr
}

real_data[1,'FDR_tmp'] = 0
real_data$FDR <- real_data$FDR_tmp

for (i in nrow(real_data):2){
    if (real_data[i,'FDR'] < real_data[i-1, 'FDR']){
        real_data[i-1,'FDR'] <- real_data[i, 'FDR']
    }
}
write.table(real_data, file = "methylation_vs_parental_bender_main_cellcount_corrected.scaled.v3.p0.05.FDR22.txt", sep = "\t", quote = F, col.names = NA)
real_data <- real_data[real_data$FDR < 0.05,]
write.table(real_data, file = "methylation_vs_parental_bender_main_cellcount_corrected.scaled.v3.FDR0.05.txt", sep = "\t", quote = F, col.names = NA)
