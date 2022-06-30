library(ieugwasr)
setwd("C:/Users/Dasha/work/UMCG/data/SV_GWAS/v4")
dsv <- read.delim("dSV_meta_all_5e-08.sorted.rsids.genes.gwas_annot.txt.clumped_0.1.txt", sep = "\t", as.is = T)
dsv <- dsv[dsv$FDR.0.05 == 1,]
dsv$proxies <- NA
g <- read.delim("food_liking.txt", sep = "\t", as.is = T)
g <- g[startsWith(g$rsID, "rs"),]
#res_proxies <- data.frame()
for (i in 1:nrow(dsv)){
  snp <- dsv[i, "rsid"]
  snp_selection <- unique(g[g$Chromosome == dsv[i,"SNPChr"],"rsID"])
  
    res <- ld_matrix(c(snp, snp_selection))
    print(i)
    
    if(!is.null(nrow(res))){
    if (nrow(res) > 1){
      proxies <- colnames(res)[which(res[1,] > 0.8)]
      proxies <- proxies[!startsWith(proxies, snp)]
      cat(snp, proxies)
    }
    }

}
