setwd("/Users/Dasha/work/UMCG/data/SV_GWAS/v4/")
f <- read.delim("clumped_vSV_res_5e-08.per_cohort_p.txt", sep = "\t", header = T, as.is = T, check.names = F)
f$Meta_replicated_in_mulitple_cohorts <- FALSE
for (i in 1:nrow(f)){
  print(i)
  print(f[i,"pvalues_per_cohort"])
  if (sum(na.omit(as.numeric(unlist(strsplit(f[i,"pvalues_per_cohort"], ",")))) < 0.05) > 1) f[i, "Meta_replicated_in_mulitple_cohorts"] <- TRUE
}

write.table(f, file = "clumped_vSV_res_5e-08.per_cohort_p.v2.txt",  sep = "\t", quote = F, col.names = NA)
