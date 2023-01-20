
setwd("/Users/Dasha/work/UMCG/data/NEXT/metabolites/")


fname <- "metabo_Baby_B.filtered.txt"
d <- as.data.frame(t(read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)))

pca <- prcomp(d, center = TRUE,scale = TRUE)
summary(pca)

# PC104 PCs explain 90% of variation
5e-08/104

# Mother P12: 4.807692e-10

fname <- "metabo_Mother_B.filtered.txt"
d <- as.data.frame(t(read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)))

pca <- prcomp(d, center = TRUE,scale = TRUE)
summary(pca)

# PC110 PCs explain 90% of variation
5e-08/110

# Mother P12: 4.545455e-10

fname <- "metabo_Mother_P12.filtered.txt"
d <- as.data.frame(t(read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)))

pca <- prcomp(d, center = TRUE,scale = TRUE)
summary(pca)

# PC104 PCs explain 90% of variation
5e-08/104

# Mother P12: 4.807692e-10

fname <- "metabo_Mother_P28.filtered.txt"
d <- as.data.frame(t(read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)))

pca <- prcomp(d, center = TRUE,scale = TRUE)
summary(pca)

# 113 PCs explain 90% of variation
5e-08/113

# Mother P28: 4.424779e-10