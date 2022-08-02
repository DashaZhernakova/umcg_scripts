library(vegan)



d <- read.delim("/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
#d <- read.delim("/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/CVD3_olinkNormal_1447_LLDsamples_ProtNames.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)

pheno_path <- paste0(config$basedir_path, "/", config$pheno_path)
pheno0 <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
clusters <- read.delim("/Users/Dasha/work/UMCG/data/gender_differences/omics/results/results/tables/explained_variance_nmr_clusters.txt",header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

d <- scale(na.omit(d))
pheno0 <- na.omit(pheno0)
ids <- intersect(row.names(d), row.names(pheno0))
for (i in 1:4){
  lipids <- row.names(clusters[clusters$cluster==i,])
  Distance = vegdist(d[ids,lipids],  method = "euclidean" )
  Formula = as.formula("Distance ~ age + gender_F1M2 + age:gender_F1M2")

  print(adonis2(formula = Formula , data = pheno0[ids,]  , na.action = na.exclude, permutations = 10000 ) )
}

lipids <- colnames(d)
Distance = vegdist(d[ids,lipids],  method = "euclidean" )
Formula = as.formula("Distance ~ age + gender_F1M2 + age:gender_F1M2")
print(adonis2(formula = Formula , data = pheno0[ids,]  , na.action = na.exclude, permutations = 10000 ) )

