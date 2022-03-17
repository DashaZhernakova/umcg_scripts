library(vegan)



d <- read.delim("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
d <- read.delim("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/CVD3_olinkNormal_1447_LLDsamples_ProtNames.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)

d <- na.omit(d)
ids <- intersect(row.names(d), row.names(pheno_m))
Distance = vegdist(d[ids,],  method = "euclidean" )
Formula = as.formula("Distance ~ age + gender_F1M2 + age:gender_F1M2")
adonis2(formula = Formula , data = pheno2[ids,]  , na.action = na.exclude, permutations = 100 ) 

