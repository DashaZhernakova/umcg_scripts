
setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/")
fitted_matrix <- read.delim("results/tables/nmr_scaled_fitted_all.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif <- read.delim("results/tables/nmr_scaled_summary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif_inters <- row.names(signif[signif$inter_p < 0.05,])

phenos <- c("S.VLDL.TG_p", "XS.VLDL.TG_p", "L.HDL.TG_p", "Gp", "XL.HDL.TG_p", "S.HDL.TG", "S.HDL.TG_p", "M.VLDL.C", "M.VLDL.CE", "S.VLDL.FC", "VLDL.C", "M.HDL.TG_p", "XS.VLDL.TG", "S.VLDL.PL", "S.VLDL.P", "S.VLDL.L", "VLDL.D", "TG.PG", "Serum.TG", "M.VLDL.PL", "M.VLDL.FC", "S.VLDL.TG", "M.VLDL.P", "M.VLDL.L", "M.VLDL.TG", "VLDL.TG", "DAG", "XXL.VLDL.C", "XXL.VLDL.CE", "XL.VLDL.PL", "XL.VLDL.TG", "XL.VLDL.P", "XL.VLDL.L", "XXL.VLDL.PL", "XXL.VLDL.TG", "XXL.VLDL.P", "XXL.VLDL.L", "XXL.VLDL.FC", "XL.VLDL.FC", "L.VLDL.TG", "L.VLDL.P", "L.VLDL.L", "L.VLDL.PL", "L.VLDL.FC", "XL.VLDL.C", "XL.VLDL.CE", "L.VLDL.C", "L.VLDL.CE")
phenos = c("L.HDL.C_p", "L.HDL.CE_p", "UnsatDeg", "XS.VLDL.FC_p", "IDL.FC_p", "IDL.PL_p", "XL.HDL.PL_p", "L.HDL.FC_p", "LA.FA", "FAw6.FA", "PUFA.FA", "M.HDL.C_p", "M.HDL.CE_p", "L.LDL.FC_p", "XS.VLDL.C_p", "XS.VLDL.CE_p")
phenos <- c("XL.HDL.C", "XL.HDL.CE", "XL.HDL.FC", "XL.HDL.P", "XL.HDL.L", "XL.HDL.PL", "HDL.D", "HDL.C", "HDL2.C", "L.HDL.PL", "L.HDL.P", "L.HDL.L", "L.HDL.FC", "L.HDL.C", "L.HDL.CE")
phenos <- c("ApoB.ApoA1", "Remnant.C", "ApoB", "S.VLDL.C", "S.VLDL.CE", "XS.VLDL.C", "XS.VLDL.CE", "XS.VLDL.FC", "XS.VLDL.PL", "XS.VLDL.P", "XS.VLDL.L", "M.LDL.PL", "S.LDL.PL", "S.LDL.FC", "M.LDL.FC", "M.LDL.P", "M.LDL.L", "S.LDL.P", "S.LDL.L", "FreeC", "Serum.C", "EstC", "IDL.FC", "L.LDL.FC", "IDL.C", "IDL.CE", "S.LDL.C", "S.LDL.CE", "M.LDL.C", "M.LDL.CE", "L.LDL.CE", "L.LDL.C", "LDL.C", "IDL.PL", "IDL.P", "IDL.L", "L.LDL.PL", "L.LDL.P", "L.LDL.L")
phenos <- c("XL.HDL.TG", "L.HDL.TG", "L.VLDL.FC_p", "M.VLDL.FC_p", "LA", "FAw6", "PUFA", "TotPG", "PC", "TotCho", "M.LDL.TG", "L.LDL.TG", "LDL.TG", "IDL.TG", "S.LDL.TG", "TotFA", "SFA", "MUFA", "M.HDL.TG", "HDL.TG")
phenos <- c("DAG.TG", "CLA", "CLA.FA", "Val", "Ile", "Leu", "Phe", "Tyr", "Glc", "Crea", "L.VLDL.TG_p", "XXL.VLDL.TG_p", "L.VLDL.PL_p", "XL.VLDL.PL_p", "XL.VLDL.TG_p", "S.HDL.PL_p", "S.HDL.FC_p", "S.VLDL.PL_p", "S.HDL.PL", "S.HDL.FC", "S.HDL.P", "S.HDL.L", "Cit", "Ala", "Lac", "Pyr", "M.VLDL.PL_p", "M.HDL.PL_p", "IDL.TG_p", "L.LDL.TG_p", "M.LDL.TG_p", "S.LDL.TG_p", "LDL.D", "L.LDL.PL_p", "M.LDL.FC_p", "S.LDL.FC_p", "M.LDL.PL_p", "S.LDL.PL_p", "XL.HDL.FC_p", "XL.HDL.C_p", "XL.HDL.CE_p", "L.HDL.PL_p")
pdf("../cvd_rfs/nmr_cluster10.pdf", height = 7, width = 10)
draw_multiple_fitted_lines(fitted_matrix[,phenos], signif_inters)
dev.off()



cvd_nmr <- c("Phe","Tyr","Lac","Pyr","bOHBut","MUFA","FAw3","DHA","PUFA.FA","MUFA.FA","UnsatDeg","Gp","ApoA1","ApoB","Serum.TG","XL.VLDL.P","XXL.VLDL.P","L.VLDL.P","M.VLDL.P","S.VLDL.P","XS.VLDL.P","IDL.P","L.LDL.P","M.LDL.P","S.LDL.P","L.HDL.P","M.HDL.P","HDL.D","VLDL.C","IDL.C","LDL.C","HDL.C","HDL2.C")
cvd_nmr_nointer <- cvd_nmr[! cvd_nmr %in% signif_inters]
cvd_nmr_inter <- cvd_nmr[cvd_nmr %in% signif_inters]
draw_multiple_fitted_lines(fitted_matrix[,cvd_nmr_inter], signif_inters)

fitted_matrix <- read.delim("results/tables/nmr_corrected_bmi_sbp_smk1_statins_all_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

pdf("../cvd_rfs/nmr_prot_rfs.pdf", height = 5, width = 10)
par(mfrow=c(2,4))
phenos <- c("HDL.D", "HDL.C", "M.HDL.P", "L.HDL.P", "HDL2.C", "ApoA1", "PUFA.FA") #no inter
draw_multiple_fitted_lines(fitted_matrix[,phenos], signif_inters, plot_title = paste(phenos,  collapse = ", "))
phenos <- c("Phe") #inverse, no inter
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,phenos]), signif_inters, plot_title = paste(phenos,  collapse = ", "))
phenos <- c("S.LDL.P", "M.LDL.P", "L.LDL.P", "IDL.P","LDL.C",  "IDL.C", "XS.VLDL.P", "ApoB","VLDL.C") 
draw_multiple_fitted_lines(fitted_matrix[,phenos], signif_inters, plot_title = paste(phenos,  collapse = ", "))
phenos <- c("XXL.VLDL.P", "XL.VLDL.P", "L.VLDL.P", "M.VLDL.P", "S.VLDL.P",  "Serum.TG", "MUFA") 
draw_multiple_fitted_lines(fitted_matrix[,phenos], signif_inters, plot_title = paste(phenos,  collapse = ", "))


phenos <- c("DHA", "FAw3", "Gp")
draw_multiple_fitted_lines(fitted_matrix[,phenos], signif_inters, plot_title = paste(phenos,  collapse = ", "))
phenos <- c("bOHBut", "Tyr", "UnsatDeg","MUFA.FA", "Gp", "Lac", "Pyr")
draw_multiple_fitted_lines(fitted_matrix[,phenos], signif_inters, plot_title = paste(phenos,  collapse = ", "))

fitted_matrix <- read.delim("results/tables/proteins_corrected_cellc_bmi_sbp_ldl_smk1_oc_statins_all_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)

phenos <- c("FABP4", "TFPI", "PAI", "TR-AP")
signif_inters <- phenos
draw_multiple_fitted_lines(fitted_matrix[,phenos], signif_inters, plot_title = paste(phenos,  collapse = ", "))
dev.off()


phenos_all <- c("XXL.VLDL.P", "XXL.VLDL.L", "XXL.VLDL.PL", "XXL.VLDL.C", "XXL.VLDL.CE", "XXL.VLDL.FC", "XXL.VLDL.TG", "XL.VLDL.P", "XL.VLDL.L", "XL.VLDL.PL", "XL.VLDL.C", "XL.VLDL.CE", "XL.VLDL.FC", "XL.VLDL.TG", "L.VLDL.P", "L.VLDL.L", "L.VLDL.PL", "L.VLDL.C", "L.VLDL.CE", "L.VLDL.FC", "L.VLDL.TG", "M.VLDL.P", "M.VLDL.L", "M.VLDL.PL", "M.VLDL.C", "M.VLDL.CE", "M.VLDL.FC", "M.VLDL.TG", "S.VLDL.P", "S.VLDL.L", "S.VLDL.PL", "S.VLDL.C", "S.VLDL.CE", "S.VLDL.FC", "S.VLDL.TG", "XS.VLDL.TG", "VLDL.C", "Remnant.C", "Serum.TG", "VLDL.TG", "DAG", "MUFA", "SFA", "CLA.FA", "XS.VLDL.P", "XS.VLDL.L", "XS.VLDL.PL", "XS.VLDL.C", "XS.VLDL.CE", "XS.VLDL.FC", "IDL.P", "IDL.L", "IDL.PL", "IDL.C", "IDL.CE", "IDL.FC", "IDL.TG", "L.LDL.P", "L.LDL.L", "L.LDL.PL", "L.LDL.C", "L.LDL.CE", "L.LDL.FC", "M.LDL.P", "M.LDL.L", "M.LDL.PL", "M.LDL.C", "M.LDL.CE", "M.LDL.FC", "S.LDL.P", "S.LDL.L", "S.LDL.PL", "S.LDL.C", "S.LDL.CE", "S.LDL.FC", "L.LDL.C_p", "L.LDL.CE_p", "M.LDL.C_p", "M.LDL.CE_p", "S.LDL.C_p", "S.LDL.CE_p", "S.HDL.C_p", "S.HDL.CE_p", "Serum.C", "LDL.C", "EstC", "FreeC", "ApoB", "ApoB.ApoA1", "TotFA", "LA", "CLA", "FAw3", "FAw6", "PUFA", "XL.HDL.P", "XL.HDL.L", "XL.HDL.PL", "XL.HDL.C", "XL.HDL.CE", "XL.HDL.FC", "XL.HDL.TG", "L.HDL.P", "L.HDL.L", "L.HDL.PL", "L.HDL.C", "L.HDL.CE", "L.HDL.FC", "L.HDL.TG", "M.HDL.P", "M.HDL.L", "M.HDL.PL", "M.HDL.C", "M.HDL.CE", "M.HDL.FC", "S.HDL.P", "S.HDL.L", "S.HDL.PL", "XS.VLDL.FC_p", "HDL.C", "HDL2.C", "DHA", "DHA.FA", "Phe", "AcAce", "Crea", "L.LDL.PL_p", "L.LDL.FC_p", "M.LDL.PL_p", "M.LDL.FC_p", "S.LDL.PL_p", "S.LDL.FC_p", "M.HDL.C_p", "M.HDL.CE_p", "S.HDL.PL_p", "TotPG", "PC", "SM", "TotCho", "ApoA1", "Glc", "Gln", "His", "Ile", "Leu", "Val", "Tyr", "bOHBut", "Alb", "L.LDL.TG", "M.LDL.TG", "S.LDL.TG", "M.HDL.TG", "S.HDL.FC", "S.HDL.TG", "XXL.VLDL.PL_p", "XXL.VLDL.C_p", "XXL.VLDL.CE_p", "XXL.VLDL.FC_p", "XXL.VLDL.TG_p", "XL.VLDL.PL_p", "XL.VLDL.C_p", "XL.VLDL.CE_p", "XL.VLDL.FC_p", "XL.VLDL.TG_p", "L.VLDL.PL_p", "L.VLDL.C_p", "L.VLDL.CE_p", "L.VLDL.FC_p", "L.VLDL.TG_p", "M.VLDL.PL_p", "S.LDL.TG_p", "S.HDL.FC_p", "S.HDL.TG_p", "HDL3.C", "LDL.TG", "HDL.TG", "TG.PG", "FALen", "UnsatDeg", "LA.FA", "FAw3.FA", "FAw6.FA", "PUFA.FA", "MUFA.FA", "SFA.FA", "Lac", "Pyr", "Cit", "Ala", "Ace", "Gp", "S.HDL.C", "S.HDL.CE", "M.VLDL.C_p", "M.VLDL.CE_p", "M.VLDL.FC_p", "S.VLDL.C_p", "S.VLDL.CE_p", "S.VLDL.FC_p", "XS.VLDL.PL_p", "M.VLDL.TG_p", "S.VLDL.PL_p", "S.VLDL.TG_p", "XS.VLDL.C_p", "XS.VLDL.CE_p", "XS.VLDL.TG_p", "IDL.PL_p", "IDL.C_p", "IDL.CE_p", "IDL.FC_p", "IDL.TG_p", "L.LDL.TG_p", "M.LDL.TG_p", "XL.HDL.FC_p", "XL.HDL.TG_p", "L.HDL.PL_p", "L.HDL.C_p", "L.HDL.CE_p", "LDL.D", "HDL.D", "XL.HDL.PL_p", "XL.HDL.C_p", "XL.HDL.CE_p", "L.HDL.FC_p", "L.HDL.TG_p", "M.HDL.PL_p", "M.HDL.FC_p", "M.HDL.TG_p", "VLDL.D", "DAG.TG")
pdf("../cvd_rfs/nmr_clusters.pdf", height = 10, width = 10)
par(mfrow=c(5,5))
for (p in phenos_all){
  draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,p]), plot_title = p)
}
dev.off()


setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/")

fitted_matrix <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_summary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
signif_inters <- row.names(signif[signif$inter_p_adj_bonferroni < 0.05,])
clusters <- read.delim("results/tables/nmr_cluster_nums.txt", header = T,  sep = "\t", as.is = T, check.names = F)
pdf("../cvd_rfs/nmr_clusters_v2.pdf", height = 10, width = 10)
par(mfrow=c(3,3))
clusters[clusters$cluster_new == "2a", 2] = 7
clusters[clusters$cluster_new == "6", 2] = 7
cnt <- 1
for (p in c("1", "2", "3", "4", "7")){
  lipids <- clusters[clusters$cluster_new == p, 1]
  if (! p == "4"){
    lipids <- lipids[lipids %in% signif_inters]
  }
  cat(p, length(lipids), "\n", sep = " ")
  draw_multiple_fitted_lines(fitted_matrix[,lipids], signif_inters, plot_title = paste0("cluster ", cnt, "\nN = ", length(lipids)))
  cnt <- cnt + 1
}
dev.off()




