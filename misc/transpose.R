args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile = args[2]
d <- read.delim(infile, header = F, row.names = 1,sep = "\t")
d2 <- as.data.frame(t(d))
write.table(d2, file = outfile, sep = "\t", quote = F, row.names = F)