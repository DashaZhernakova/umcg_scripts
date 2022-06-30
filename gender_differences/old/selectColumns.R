args <- commandArgs(trailingOnly = TRUE)


d <- read.delim(args[1], header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)

cols <- read.delim(args[2], header = F, sep = "\t", as.is = T, check.names = F)
cols <- cols[,1]

d2 <- d[,cols]
write.table(d2, file = args[3], sep = "\t", quote = F, col.names = NA)

