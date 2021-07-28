library(CVrisk)

d <- read.delim("age_gender_bmi_glu_smk_lipids_sbp_t2d_antihyper_all_LL_with_children.txt", check.names = F, as.is = T, sep = "\t", header = T, row.names = 1)
d <- d[d$age < 74 & d$age > 30,]
d2 <- sapply(d, function(x) as.numeric(as.character(x)))
row.names(d2) <- row.names(d)

sex <- rep(NA, nrow(d2))
sex[d2[,"gender_F1M2"] == 1] = "female"
sex[d2[,"gender_F1M2"] == 2] = "male"
d2 <- subset(d2, select = c(-GLU, -gender_F1M2))

frs <- matrix(nrow = nrow(d2), ncol = 2)
colnames(frs) = c("sample", "FRS")
frs[,1] <- row.names(d2)
for (i in 1:nrow(d2)){
    frs [i,2] <- ascvd_10y_frs(gender = sex[i], age = d2[i,1], hdl = d2[i, 5], totchol = d2[i, 4], sbp = d2[i, 7],  bp_med = d2[i,9], smoker = d2[i, 3], diabetes = d2[i,8])
}
write.table(frs, file = paste0(out_table_path, "_summary.txt"), sep = "\t", quote = F, col.names = NA)
