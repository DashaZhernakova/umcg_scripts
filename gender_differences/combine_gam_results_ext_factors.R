library(dplyr)
library(corrplot)
phenos <- c("AF", "ALT", "AST", "CA", "NAA", "FOS", "GLU", "DBP", "CHO", "HDC", "LDC","TGL")

ptable <- data.frame()
stable <- data.frame()
pval_table <- data.frame()

cnt = 1
for (ph in phenos){
    print(ph)
    s <- read.delim(paste0(ph,".s.table.txt"), sep = "\t", as.is = T, check.names = F, header = T)
    p <- read.delim(paste0(ph,".p.table.txt"), sep = "\t", as.is = T, check.names = F, header = T)
    
    colnames(s) <- c("V1", paste(colnames(s)[2:ncol(s)], rep(ph, 4), sep = "_"))
    colnames(p) <- c("V1", paste(colnames(p)[2:ncol(p)], rep(ph, 4), sep = "_"))
    if (cnt == 1){
        stable <- s
        ptable <- p
        #pval_table <- cbind(c(s[,1], p[,1]),c(s[,5], p[,5]))
        #colnames(pval_table) <- c("V1", ph)
    } else {
        stable <- full_join(stable, s, by="V1")
        ptable <- full_join(ptable, p, by = "V1")
        #pval_table <- full_join(pval_table, rbind(s[,c(1,5)], p[,c(1,5)]), by = "V1")
    }
    cnt = cnt +1
}

write.table(stable, file = "stable.txt", quote = F, col.names = NA, sep = "\t")

ppval_table <- ptable[,seq(5,ncol(ptable), 4)]
spval_table <- stable[,seq(5,ncol(stable), 4)]
colnames(ppval_table) <- phenos
colnames(spval_table) <- phenos
res <- rbind(ppval_table, spval_table)
row.names(res) <- c(ptable$V1, stable$V1)


write.table(res, file = "pval_table.txt", quote = F, col.names = NA, sep = "\t")
pval_bins <- as.data.frame(apply(res, 2, cut, c(-Inf, 0.0001, 0.05, 1), labels=c(1,0.5,0)))
row.names(pval_bins) <- row.names(res)
write.table(pval_bins, file = "pval_table_bins.txt", quote = F, col.names = NA, sep = "\t")

pval_bins <- sapply(pval_bins, function(x) as.numeric(as.character(x)))
row.names(pval_bins) <-row.names(res)
pdf("pval_bins.corrplot.pdf", height = 10, width = 10)
corrplot(pval_bins)
dev.off()

#TODO: rename rows!