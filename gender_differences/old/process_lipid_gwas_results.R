
chr= 22
gender = "women"
age_groups = c("below_40", "above_60")
column <- "BETA"
column <- "P"
column <- c("BETA", "P")
#  for (pheno in c("HDC", "LDC", "RCHO", "TGL")){
pheno="LDC"
    all_chr <- data.frame()
    for (chr in seq(1,22)){
      print(chr)
      gender = "women"
      inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[1], pheno)
      d1 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      
      d1[d1$ID == ".", "ID"] <- paste0(d1[d1$ID == ".","#CHROM"], ":", d1[d1$ID == ".","POS"])
      d1 <- d1[!duplicated(d1$ID),]
      rownames(d1) <- d1$ID
      
      inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[2], pheno)
      d2 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      
      d2[d2$ID == ".", "ID"] <- paste0(d2[d2$ID == ".","#CHROM"], ":" , d2[d2$ID == ".","POS"])
      d2 <- d2[!duplicated(d2$ID),]
      rownames(d2) <- d2$ID
      
      #gender = "men"
      #inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[1], pheno)
      #d3 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      #d3 <- d3[!duplicated(d3$POS),]
      #rownames(d3) <- d3$POS
      
      #inpath <- sprintf("chr%s.assoc.%s_%s_unrel.%s.glm.linear.gz", chr, gender, age_groups[2], pheno)
      #d4 <- read.delim(gzfile(inpath), header = T, as.is = T, check.names = F)
      #d4 <- d4[!duplicated(d4$POS),]
      #rownames(d4) <- d4$POS
      
      #all_snps <- intersect(intersect(row.names(d1), row.names(d2)), intersect(row.names(d3), row.names(d4)))
      all_snps <- intersect(row.names(d1), row.names(d2))
      betas <- cbind(d1[all_snps,column[1]], d2[all_snps,column[1]], d1[all_snps,column[2]], d2[all_snps,column[2]])
      row.names(betas) <- all_snps
      colnames(betas) <- c( paste0("women_below_40_", column[1]), paste0("women_above_60_", column[1]), 
                            paste0("women_below_40_", column[2]), paste0("women_above_60_", column[2]) )
      all_chr <- rbind(all_chr, betas)
    }
    write.table(all_chr, file = paste0("all_chr.", pheno, "corrected.merged.txt"), sep = "\t", quote = F, col.names = NA)
  #}
#}

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}


setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/all_LL/lipid_gwas/")
d <- read.delim(gzfile("all_chr.LDC.corrected.merged.1e-4.txt"), header = T, as.is = T, check.names = F)

pdf("ldc_betas_plot.pdf")
#plot(d$women_below_40_BETA, d$women_above_60_BETA)
ggplot(d, aes(x = women_below_40_BETA, y = women_above_60_BETA)) + 
  geom_point() + xlim(-0.5,0.5) + ylim(-0.5, 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "grey") + geom_abline(intercept = 0, slope = -1, color = "grey")
dev.off()
