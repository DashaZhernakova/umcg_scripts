
chr= 22
gender = "women"
age_groups = c("below_40", "above_60")
pheno = "LDC"
column <- "BETA"
#column <- "P"

clumped_snps <- read.delim(sprintf("%s_%s_clumped_all.txt", pheno, gender), as.is = T, check.names = F)
clumped_snps <- clumped_snps[,1]

for (column in c("BETA", "p.value")){
  #for (pheno in c("LDC")){
    all_chr <- data.frame()
    for (chr in seq(1,22)){
      print(chr)
      gender = "women"
      inpath <- sprintf("%s_%s_%s_gwas/p-values.chr.%s", pheno, gender, age_groups[1], chr )
      d1 <- read.delim(inpath, header = T, as.is = T, check.names = F, sep = " ")
      d1 <- d1[!duplicated(d1$POS),]
      rownames(d1) <- d1$SNPID
      
      inpath <- sprintf("%s_%s_%s_gwas/p-values.chr.%s", pheno, gender, age_groups[2], chr )
      d2 <- read.delim(inpath, header = T, as.is = T, check.names = F, sep = " ")
      d2 <- d2[!duplicated(d2$POS),]
      rownames(d2) <- d2$SNPID
      
      # gender = "men"
      # inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[1], pheno)
      # d3 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      # d3 <- d3[!duplicated(d3$POS),]
      # rownames(d3) <- d3$POS
      # 
      # inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[2], pheno)
      # d4 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      # d4 <- d4[!duplicated(d4$POS),]
      # rownames(d4) <- d4$POS
      
      #all_snps <- intersect(row.names(d1), row.names(d2))
      all_snps <- clumped_snps[grepl(paste0("^",chr, ":"), clumped_snps)]
      
      betas <- cbind(d1[all_snps,column], d2[all_snps,column])
      row.names(betas) <- paste0(chr, ":", all_snps)
      colnames(betas) <- c("women_below_40", "women_above_60")
      
      all_chr <- rbind(all_chr, betas)
    }
    write.table(all_chr, file = paste0("all_chr.", pheno, ".merged_", column, ".clumped.txt"), sep = "\t", quote = F, col.names = NA)
  #}
}

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}


setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/all_LL/lipid_gwas/")
d <- read.delim(gzfile("all_chr.LDC.merged_p.value.clumped.txt.gz"), header = T, as.is = T, check.names = F)

png("ldc_clumped_betas_plot.png", units = 'in', res = 400)
plot(d$women_below_40, d$women_above_60)
dev.off()