
chr= 22
gender = "women"
age_groups = c("below_40", "above_60")
pheno = "CHO"
column <- "BETA"
column <- "P"
for (column in c("BETA", "P")){
  for (pheno in c("HDC", "LDC", "RCHO", "TGL")){
    all_chr <- data.frame()
    for (chr in seq(1,22)){
      print(chr)
      gender = "women"
      inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[1], pheno)
      d1 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      d1 <- d1[!duplicated(d1$POS),]
      rownames(d1) <- d1$POS
      
      inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[2], pheno)
      d2 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      d2 <- d2[!duplicated(d2$POS),]
      rownames(d2) <- d2$POS
      
      gender = "men"
      inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[1], pheno)
      d3 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      d3 <- d3[!duplicated(d3$POS),]
      rownames(d3) <- d3$POS
      
      inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[2], pheno)
      d4 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      d4 <- d4[!duplicated(d4$POS),]
      rownames(d4) <- d4$POS
      
      all_snps <- intersect(intersect(row.names(d1), row.names(d2)), intersect(row.names(d3), row.names(d4)))
      
      betas <- cbind(d1[all_snps,column], d2[all_snps,column], d3[all_snps,column], d4[all_snps,column])
      row.names(betas) <- paste0(chr, ":", all_snps)
      colnames(betas) <- c("women_below_40", "women_above_60", "men_below_40", "men_above_60")
      
      all_chr <- rbind(all_chr, betas)
    }
    write.table(all_chr, file = paste0("all_chr.", pheno, ".merged_", column, ".txt"), sep = "\t", quote = F, col.names = NA)
  }
}

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}


setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/all_LL/lipid_gwas/")
d <- read.delim(gzfile("all_chr.CHO.merged_P.txt.gz"), header = T, as.is = T, check.names = F)

png("cho_betas_plot.png", units = 'in', res = 400)
plot(d$women_below_40, d$women_above_60)
dev.off()