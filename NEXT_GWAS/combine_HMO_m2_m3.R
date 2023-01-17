m2 <- read.delim("../HMO_GWAS/HMO_M2.txt", header = T, row.names = 1, sep = "\t", as.is = T)
m3 <- read.delim("../HMO_GWAS/HMO_M3.txt", header = T, row.names = 1, sep = "\t", as.is = T)

overlap <- intersect(row.names(m2), row.names(m3))

m2o <- m2[overlap,]
m3o <- m3[overlap,]

m23 <- (m2o + m3o)/2
diag(cor(m2o, m3o))
