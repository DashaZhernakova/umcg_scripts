library(corrplot)
phenos <- c("AF", "ALT", "AST","CA", "FOS", "NAA", "GLU",  "DBP", "CHO", "HDC", "LDC", "TGL","GR")
path <- "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/factors/results/"

bootstrap_res <- data.frame()
n_boot <- 10

cnt = 1
for (p in phenos){
  p_table <- read.delim(paste0(path, p, ".p.table.txt"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
  s_table <- read.delim(paste0(path, p, ".s.table.txt"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

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
}

pdf("lifestyle_factors_pval_corrplot.pdf")
corrplot(as.matrix(1-res),is.corr = F, tl.col = "black")


pdf("lifestyle_factors_pvalbins_corrplot.pdf")
cuts <- apply(res, 2, function(x) {4-as.numeric(cut(x, c(-Inf,0.001, 0.05,1), labels=1:3))})

row.names(cuts) <- row.names(res)
corrplot(cuts,is.corr = F, tl.col = "black")
dev.off()



############
bsum <- data.frame()
n_boot <- 10
#phenos <- c("AF", "ALT", "AST","CA", "FOS", "NAA", "GLU",  "DBP", "CHO", "HDC",  "TGL")

for (b in 1:n_boot){
  cnt = 1
  for (p in phenos){
    p_table <- read.delim(paste0(path, p, "_bootstrap_", b, ".p.table.txt"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
    s_table <- read.delim(paste0(path, p, "_bootstrap_",b , ".s.table.txt"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
    
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
  }
  cuts <- apply(res, 2, function(x) {ifelse(x < 0.05, 1, 0)})
  cuts[is.na(cuts)] <- 0
  if (b == 1){
    bsum <- cuts
  } else {
    bsum <- bsum + cuts
  }
}

pdf("lifestyle_factors_bootstrap_corrplot.pdf")
corrplot(bsum,is.corr = F, tl.col = "black")
dev.off()

bsum2 <- t(bsum)
bsum2[bsum2 > 1 & bsum2 < 9] <- NA
cols <- c("Intercept", "sex", "smoking", "smoking:sex1.1", "smoking:sex0.2", "smoking:sex1.2", "med", "med:sex1.1", "med:sex0.2", "med:sex1.2", "age", "age:sex", "stress_1y", "stress_chronic", "phys_activity_total", "phys_activity_intensive", "diet", "alcohol", "stress_1y:age:sex", "stress_chronic:age:sex", "phys_activity_total:age:sex", "phys_activity_intensive:age:sex", "diet:age:sex", "alcohol:age:sex", "stress_1y:age", "stress_chronic:age", "phys_activity_total:age", "phys_activity_intensive:age", "diet:age", "alcohol:age", "stress_1y:sex", "stress_chronic:sex", "phys_activity_total:sex", "phys_activity_intensive:sex", "diet:sex", "alcohol:sex", "age:smoking", "age:smoking:sex0.1", "age:smoking:sex1.1", "age:smoking:sex0.2", "age:smoking:sex1.2", "age:med", "age:med:sex0.1", "age:med:sex1.1", "age:med:sex0.2", "age:med:sex1.2")
colnames(bsum2) <- cols
pdf("lifestyle_factors_bootstrap_0.05_filtered_corrplot.pdf")
corrplot(bsum2,is.corr = F, tl.col = "black", tl.cex = 0.8)
dev.off()
