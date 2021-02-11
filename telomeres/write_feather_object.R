library(feather)
setwd('/groups/umcg-lld/scr01/dasha/methylation/permutation_FDR')

traits_path <- "../ILLUMINA450K_All2Blood_mValue.Selection.Dasen.QuantileNormalized.LLDsubset.txt.gz"
gte_path <- "../methylation_gte.txt"
traits0 <- as.data.frame(t(read.table(gzfile(traits_path), header = F, row.names = 1, sep = "\t", as.is = T, check.names = F)))
gte <- read.table(gte_path, sep = "\t", as.is = T, check.names = F)
gte_m <- gte[match(traits0[,1], gte[, 2], nomatch = 0 ),]
traits <- traits0[match(gte_m[, 2], traits0[,1], nomatch = 0), ]
all(traits[,1] == gte_m[, 2])
traits[,1] <- gte_m[, 1]
colnames(traits)[1] <- "sample"
 traits2 <- data.frame(apply(traits, 2, function(x) as.numeric(as.character(x)))
 traits2$sample <- traits$sample
 traits2[, -c(1)] <- scale(traits2[, -c(1)])
write_feather(traits2, "../ILLUMINA450K_All2Blood_mValue.Selection.Dasen.QuantileNormalized.LLDsubset.LLD_ids.scaled.feather")
