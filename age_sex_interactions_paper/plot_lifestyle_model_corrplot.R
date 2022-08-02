library(corrplot)
phenos <- c("AF", "ALT", "AST","CA", "FOS", "NAA",  "SBP", "CHO",  "LDC", "TGL", "HALB","GR", "ER", "HT", "BALB", "UZ")
#path <- "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/factors/results/"
setwd("/Users/Dasha/work/UMCG/data/gender_differences/omics/all_LL/lifestyle_factors/")
path <- "tables/"

bootstrap_res <- data.frame()
n_boot <- 50


############
bsum <- data.frame()
n_boot <- 50
n_boot_per_pheno <- matrix(data = 0, nrow = length(phenos), ncol = 1)
row.names(n_boot_per_pheno) <- phenos
for (b in 1:n_boot){
  cnt = 1
  for (p in phenos){
    tryCatch({
    p_table <- read.delim(paste0(path, p, "_bootstrap_", b, ".p.table.txt"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
    s_table <- read.delim(paste0(path, p, "_bootstrap_",b , ".s.table.txt"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
    
    n_boot_per_pheno[p,1] <- n_boot_per_pheno[p,1] + 1
    
    rows <- c(row.names(p_table), row.names(s_table))
    res_line <- c(p_table[,4], s_table[,4])
    if (cnt == 1){
      
      res <- data.frame(matrix(nrow = length(rows), ncol = length(phenos)))
      row.names(res) <- rows
      colnames(res) <- phenos
    }
    
    m <- match(row.names(res), rows, nomatch = 0)
    res[match(rows,row.names(res), nomatch = 0),p] <- res_line[match(row.names(res), rows, nomatch = 0)]
    
    cnt <- cnt +1
    
  },error=function(e) {print(paste0("ERROR! ", p, "_bootstrap_", b, ".p.table.txt"))})
  }
  cuts <- apply(res, 2, function(x) {ifelse(x < 0.05, 1, 0)})
  cuts[is.na(cuts)] <- 0
  if (b == 1){
    bsum <- cuts
  } else {
    bsum <- bsum + cuts
  }
}


# 
# bfrac <- apply(bsum, 1, function(x) x/n_boot_per_pheno)
# row.names(bfrac) <- phenos
# pdf("lifestyle_factors_bootstrap50_corrplot3.pdf")
# corrplot(bsum,is.corr = F, tl.col = "black", tl.cex = 0.7, col.lim = c(0,50), 
#          col=colorRampPalette(c("#FFFFFF","#08519C"))(50))
#          
# dev.off()

factors <- c("diet", "alcohol", "phys_activity_total"  ,"stress_chronic", "sleep_duration")
pdf("lifestyle_factors_bootstrap30_corrplot_split.discr3.pdf", useDingbats = F)
par(mfrow=c(6,1))
subs <- bsum[c("gender_F1M22", "s(age)", "s(age):gender_F1M21", "s(age):gender_F1M22"),]
corrplot(subs,is.corr = F, tl.col = "black", tl.cex = 0.7, col.lim = c(0,30), 
         col=colorRampPalette(c("#FFFFFF","#08519C"))(30))
for (f in factors){
  print(f)
  ord <- c(paste0("s(", f, ")"), paste0("ti(", f, ",age)"), paste0("s(", f, "):gender_F1M21"), paste0("s(", f, "):gender_F1M22"), paste0("ti(", f, ",age):gender_F1M21"), paste0("ti(", f, ",age):gender_F1M22"))
  subs <- bsum[grep(f,row.names(bsum)),]
  corrplot(subs[ord,],is.corr = F, tl.col = "black", tl.cex = 0.7, col.lim = c(0,30), 
           col=colorRampPalette(c("#FFFFFF","#08519C"))(30))
  
}

bin_fac <- c("smoking", "med")
for (f in factors){
  print(f)
  ord <- c(paste0(f, "1"), paste0("s(age):", f, "0"), paste0("s(age):", f, "1"), paste0("s(age):interaction(", f, ", gender_F1M2)0.1") , paste0("s(age):interaction(", f, ", gender_F1M2)1.1") , paste0("s(age):interaction(", f, ", gender_F1M2)0.2") , paste0("s(age):interaction(", f, ", gender_F1M2)1.2"))
  subs <- bsum[grep(f,row.names(bsum)),]
  corrplot(subs,is.corr = F, tl.col = "black", tl.cex = 0.7, col.lim = c(0,30), 
           col=colorRampPalette(c("#FFFFFF","#08519C"))(30))
  
}
dev.off()

