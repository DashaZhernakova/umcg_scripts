args <- commandArgs(trailingOnly = TRUE)
fam_f <- args[1]
perm_fam_f <- args[2]

fam <- read.table(fam_f, header = F, as.is = T, check.names = F)
fam_perm <- fam[sample(nrow(fam)),]
write.table(fam_perm, file = perm_fam_f, sep = "\t", quote = F, col.names = NA)