library(vegan)



d <- read.delim("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
d <- read.delim("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/CVD3_olinkNormal_1447_LLDsamples_ProtNames.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)

pheno_path <- paste0(config$basedir_path, "/", config$pheno_path)
pheno0 <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

d <- scale(na.omit(d))
pheno0 <- na.omit(pheno0)
ids <- intersect(row.names(d), row.names(pheno0))
Distance = vegdist(d[ids,],  method = "euclidean" )
Formula = as.formula("Distance ~ age + gender_F1M2 + age:gender_F1M2 + I(age^2) + I(age^3) + I(age^2):gender_F1M2 + I(age^3):gender_F1M2")
adonis2(formula = Formula , data = pheno0[ids,]  , na.action = na.exclude, permutations = 100 ) 

