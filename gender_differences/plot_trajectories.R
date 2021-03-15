
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
pdf("../cvd_rfs/nmr_cluster8.pdf", height = 7, width = 10)
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
