args <- commandArgs(trailingOnly = TRUE)
d <- read.delim(args[1], header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)
d_norm <- apply(d, 2, function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x))))
write.table(d_norm, file = args[2],  sep = "\t", quote = F, col.names = NA)