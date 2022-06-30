setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/")

  # d1 <- read.delim("data/CVD3_olinkNormal_1447_LLDsamples_ProtNames.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
  # d2 <- read.delim("data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
  # m <- merge_match(d1,d2)
  # d3 <- as.data.frame(pheno[,"smk1"])
  # rownames(d3) <- rownames(pheno)
  # colnames(d3) <- "smk1"
  # m2 <- merge_match(m,d3)
  # colnames(m2)[ncol(m2)] = "smk1"
  # write.table(m2, file = "data/cvd_risk_factors_33nmr.txt", sep = "\t", quote = F, col.names = NA)
  # 
  # cvd_nmr <- c("Phe","Tyr","Lac","Pyr","bOHBut","MUFA","FAw3","DHA","PUFA.FA","MUFA.FA","UnsatDeg","Gp","ApoA1","ApoB","Serum.TG","XL.VLDL.P","XXL.VLDL.P","L.VLDL.P","M.VLDL.P","S.VLDL.P","XS.VLDL.P","IDL.P","L.LDL.P","M.LDL.P","S.LDL.P","L.HDL.P","M.HDL.P","HDL.D","VLDL.C","IDL.C","LDL.C","HDL.C","HDL2.C")
  # d2 <- d2[,colnames(d2) %in% cvd_nmr]


d <- read.delim("data/cvd_risk_factors_33nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
d <- na.omit(d)
cormat <- cor(d)
#hc <- hclust(as.dist(1-cormat))
#clusters <- cutree(hc, h=0.89)

D1 <- as.dist(1-cormat)
C1 <- hclust(as.dist(1-cormat))
pdf("../cvd_rfs/heatmap_test.pdf", width = 10, height = 10)
h <- pheatmap(D1, cluster_rows = C1, cluster_cols = C1, labels_row = C1$labels[C1$order], labels_col = C1$labels[C1$order], clustering_method = 'complete', clustering_distance_rows = D1, clustering_distance_cols = D1, fontsize = 8)
dev.off()


inter_subs <- c("Gln","S.VLDL.CE","Tyr","ApoB.ApoA1","Remnant.C","S.VLDL.C","IDL.CE","IDL.C","S.LDL.PL","ApoB","L.LDL.CE_p","XS.VLDL.C","XS.VLDL.CE","IDL.CE_p","EstC","Serum.C","XS.VLDL.L","M.LDL.PL","XS.VLDL.P","M.VLDL.CE","FreeC","IDL.C_p","S.HDL.PL_p","VLDL.C","XXL.VLDL.CE","TotFA","S.LDL.FC","L.LDL.PL","S.LDL.P","S.LDL.L","L.LDL.CE","L.LDL.PL_p","L.LDL.P","IDL.L","L.LDL.C","M.LDL.L","M.LDL.P","L.LDL.L","XS.VLDL.PL","L.LDL.FC_p","M.LDL.C","M.LDL.CE","IDL.P","S.VLDL.FC","LDL.C","XS.VLDL.FC","M.LDL.FC","M.LDL.TG_p","S.LDL.C","S.HDL.CE","M.VLDL.C","XXL.VLDL.C","S.VLDL.L","XL.HDL.PL_p","M.LDL.FC_p","MUFA","M.LDL.CE_p","IDL.PL","S.LDL.CE","Alb","S.LDL.FC_p","L.VLDL.CE","S.HDL.CE_p","S.HDL.C","IDL.FC","L.VLDL.C","IDL.PL_p","L.LDL.C_p","PUFA","S.VLDL.PL_p","S.VLDL.P","S.VLDL.PL","SM","XL.HDL.C_p","Val","FAw6","L.LDL.FC","SFA","XXL.VLDL.PL","Pyr","M.LDL.PL_p","M.VLDL.FC","L.VLDL.FC","XL.VLDL.CE","L.VLDL.FC_p","XL.VLDL.C","XXL.VLDL.L","XXL.VLDL.P","L.LDL.TG_p","XXL.VLDL.TG","S.LDL.CE_p","XL.HDL.CE_p","M.HDL.CE_p","XL.VLDL.FC","L.VLDL.L","M.LDL.C_p","S.LDL.PL_p","S.HDL.C_p","M.VLDL.PL","L.VLDL.P","M.VLDL.L","XXL.VLDL.CE_p","XXL.VLDL.FC","M.VLDL.P","L.VLDL.TG","FAw3","L.VLDL.PL","M.LDL.TG","L.LDL.TG","L.VLDL.TG_p","LA","VLDL.TG","LDL.TG","L.VLDL.C_p","Glc","Leu","DAG","M.HDL.TG_p","S.LDL.C_p","CLA","M.VLDL.FC_p","IDL.TG_p","IDL.TG","Serum.TG","S.LDL.TG","His","XS.VLDL.TG","M.HDL.TG","M.VLDL.TG","XL.HDL.FC_p","XL.VLDL.L","XL.VLDL.PL","XL.HDL.TG_p","TotCho","XL.VLDL.P","HDL3.C","TotPG","LDL.D","bOHBut","M.HDL.C_p","XL.VLDL.CE_p","XL.VLDL.PL_p","XL.VLDL.TG_p","XL.VLDL.TG","L.HDL.CE_p","S.LDL.TG_p","XL.VLDL.FC_p","CLA.FA","L.HDL.TG_p","DAG.TG","L.HDL.C_p","Ile","S.VLDL.TG","DHA","SFA.FA","DHA.FA","PC","HDL.TG","M.VLDL.TG_p","UnsatDeg","FAw3.FA","XL.VLDL.C_p","XXL.VLDL.C_p","M.HDL.PL_p","L.HDL.FC_p","Lac","M.VLDL.C_p","XS.VLDL.PL_p","Gp","S.HDL.TG","MUFA.FA","XXL.VLDL.PL_p","LDL-receptor","GDF-15","COL1A1","CHIT1","PAI","FABP4","TFPI","TFF3","PON3","ITGB2","TR-AP","uPA","IGFBP-2","AXL","CDH5","SCGB3A2","FAS","MMP-3","DLK-1","Ep-CAM","RETN","EGFR","CPB1","CD163","EPHB4","MB","AP-N")
d_i <- d[,colnames(d) %in% inter_subs]
cormat_i <- cor(d_i)
D1 <- as.dist(1-cormat_i)
C1 <- hclust(as.dist(1-cormat_i))
pdf("../cvd_rfs/heatmap_test_inter_subset.pdf", width = 10, height = 10)
h <- pheatmap(D1, cluster_rows = C1, cluster_cols = C1, labels_row = C1$labels[C1$order], labels_col = C1$labels[C1$order], clustering_method = 'complete', clustering_distance_rows = D1, clustering_distance_cols = D1, fontsize = 8)
dev.off()

pdf("../cvd_rfs/corrplot_test_inter_subset.pdf", width = 10, height = 10)
corrplot(cormat, method = "square", order = 'hclust', addrect = T, tl.col = "black", tl.cex = 0.5)
dev.off()

cormat2 <- cormat
cormat2[cormat2 < 0.5] <- 0
pdf("../cvd_rfs/corrplot_test_inter_subset2.pdf", width = 10, height = 10)
corrplot(cormat2, method = "square", order = 'hclust', addrect = T, tl.col = "black", tl.cex = 0.5)
dev.off()




d2 <- read.delim("data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
d2 <- na.omit(d2)
cormat_nmr <- cor(d2)
#cormat_nmr[cormat_nmr < 0.9] <- 0
D1 <- as.dist(1-cormat_nmr)
C1 <- hclust(as.dist(1-cormat_nmr))
clusters <- cutree(C1, h=1.2)
h <- pheatmap(D1, cluster_rows = C1, cluster_cols = C1, labels_row = C1$labels[C1$order], labels_col = C1$labels[C1$order], clustering_method = 'complete', clustering_distance_rows = D1, clustering_distance_cols = D1, fontsize = 8)

pdf("../cvd_rfs/corrplot_nmr.pdf", width = 10, height = 10)
corrplot(cormat_nmr, method = "square", order = 'hclust',  tl.col = "black", tl.cex = 0.2, addgrid.col = NA, cl.pos = 'n')
dev.off()



