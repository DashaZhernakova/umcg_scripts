library(reshape2)
library(dendextend)

indices = 1:ncol(traits_m)

res_rsq_table <- data.frame(matrix(nrow = length(indices), ncol = 4))
res_rsq_table_perc <- data.frame(matrix(nrow = length(indices), ncol = 4))
row.names(res_rsq_table) <- colnames(traits_m)[indices]
colnames(res_rsq_table) <- c("null_model", "sex", "age", "inter")
row.names(res_rsq_table_perc) <- colnames(traits_m)[indices]
colnames(res_rsq_table_perc) <- c("null_model", "sex", "age", "inter")
cnt = 1
cat("\nStarting the analyses\n")
for (idx in indices){
  
  trait_name <- ifelse(is.null(pheno_table), colnames(traits_m)[idx], pheno_table[pheno_table[,1] == colnames(traits_m)[idx], 2])
  cat(idx, " : ", trait_name, "\n")
  
  log_transform = FALSE
  if (colnames(traits_m)[idx] %in% pheno_to_log) log_transform = TRUE
  cat("\tLog tranform: ", log_transform, "\n")
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = "IQR", log_tr = log_transform, scale_tr = scale_transform)

  covariateslinear2 <- covariateslinear[covariateslinear != colnames(traits_m)[idx]]
  covariatesnonlinear2 <- covariatesnonlinear[covariatesnonlinear != colnames(traits_m)[idx]]
  
  exp_var_dif <- calc_explained_variance_CV(merged_tab, trait_name, covariates_linear = covariateslinear2, covariates_nonlinear = covariatesnonlinear2)
  #exp_var_dif <- c(exp_var[1], exp_var[2:4] - exp_var[1:3])
  exp_var_dif[exp_var_dif < 0] <- 0
  res_rsq_table[cnt,] <- exp_var_dif
  res_rsq_table_perc[cnt,] <- exp_var_dif/sum(exp_var_dif)
  cnt <- cnt + 1
}
write.table(res_rsq_table, file = paste0(out_basepath, "/tables/explained_variance_olink.txt"), sep = "\t", quote = F, row.names = F)
res_rsq_table <- read.delim(paste0(out_basepath, "/tables/explained_variance_nmr.txt"), sep = "\t",  as.is = T, header = 1)
row.names(res_rsq_table) <- colnames(traits_m)
# from https://stackoverflow.com/questions/67504090/add-percentage-labels-inside-bars-in-circos-barplot-in-circlize
library(circlize)
 

# NMR
setwd("/Users/Dasha/work/UMCG/data/gender_differences/omics/results/")

max_h <- max(rowSums(res_rsq_table))
t <- (res_rsq_table/max_h)
t$rest <- 1-rowSums(t)
barcolor = c("#BDBDBD", "#66C2A5", "#FC8D62", "#8DA0CB",  "#FFFFFF")
labelcolor <- c(brewer.pal(n = 8, name = "Dark2"), brewer.pal(n = 12, name = "Paired"))
#hc=reorder(hclust(dist(t)),-as.matrix(t)%*%seq(ncol(t))^2)
d <- read.delim("/Users/Dasha/work/UMCG/data/gender_differences/omics/results/results/tables/nmr_corrected_bmi_smk1_statins_bonferroni.gam_coefficients.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
fitted_matrix <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_fitted.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
#d <- read.delim("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/CVD3_olinkNormal_1447_LLDsamples_ProtNames.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
#cormat <- cor(t(d), method = "spearman", use="pairwise.complete.obs")

#hc <- hclust(as.dist(1-cormat))
hc <- hclust(dist(t(fitted_matrix)), method = "average")
#hc0 <- hc
#hc <- hc$hclust

t["axis",] <- rep(0, ncol(t))

labels=c(hc$labels[hc$order], "axis")
ord <- c(hc$order, nrow(t))
#cut=cutree(hc, k = num_k)
dend <- as.dendrogram(hc)
#dend=color_branches(as.dendrogram(hc),k=length(unique(cut)),
#                    col=labelcolor[unique(cut[labels])])
circos.clear()
pdf(paste0(out_basepath, "plots/circos_barplot_nmr.gam_hclust_eucl_fitted_cv3.pdf"),width = 20, height = 20, useDingbats = F)
#png(paste0(out_basepath, "plots/circos_barplot_olink.abs.hclust_orig_cv.png"),width = 2000, height = 2000, res = 300)

circos.par(cell.padding=c(0,0,0,0), start.degree = 0)
circos.initialize("a",xlim=c(0,nrow(t)))

circos.track(ylim=c(0,1),track.height=.2,track.margin=c(0,0),bg.border=NA,
             panel.fun=function(x,y)for(i in 1:nrow(t))circos.text(i-.5,0,labels[i],adj=c(0,.5),
                                                                   facing="clockwise",niceFacing=T,cex=1.3,col= "black"))

circos.track(ylim=c(0,1),track.height=.3,track.margin=c(0,.01),bg.border=NA,
             panel.fun=function(x,y)circos.barplot(as.matrix(t)[ord,],-.5+1:nrow(t),
                                                   col=barcolor,bar_width=1,lwd=.3,border="lightgray"))

circos.track(ylim=c(0,attr(dend,"height")),track.height=.4,track.margin=c(0,0),
             bg.border=NA,panel.fun=function(x,y)circos.dendrogram(dend))

circos.clear()
dev.off()


#
# Plot selected
#
signif <- read.delim("results/tables/nmr_corrected_bmi_smk1_statins_bonferroni_summary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
ids <- pheno_table[match(row.names(signif), pheno_table[,2], nomatch = 0),1]
row.names(signif) <- ids

signif_inters <- row.names(signif[signif$inter_p_adj_bonferroni < 0.05,])

pdf(paste0("results/plots/nmr_trajectories_selected.v2022.pdf"), width = 20, height = 20)
par(mfrow=c(5,4))

all_clusters_fitted <- data.frame()
st <- "Crea"
end <- "Leu"
lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))

st <- "S.HDL.PL_p"
end <- "S.LDL.PL_p"

lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- t(fitted_matrix[,lipids])

st <- "XL.HDL.FC"
end <- "L.HDL.CE"

lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

st <- "S.HDL.C_p"
end <- "IDL.FC"
lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

st <- "XL.HDL.C"
end <- "TotCho"

lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

st <- "L.VLDL.TG"
end <- "M.HDL.TG_p"
lipids <- labels[which(labels == st) : which(labels == end)]
draw_multiple_fitted_lines(as.data.frame(fitted_matrix[,lipids]), signif_inters, plot_title = paste0(st, " - ", end))
all_clusters_fitted <- rbind(all_clusters_fitted, t(fitted_matrix[,lipids]))

dev.off()

write.table(all_clusters_fitted, file =  "results/tables/nmr_trajectories_selected.v2022.txt", sep = "\t", quote = F, col.names = NA)
## Percentage
t <- res_rsq_table_perc
#hc=reorder(hclust(dist(t)),-as.matrix(t)%*%seq(ncol(t))^2)

barcolor = c("#BDBDBD", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3")

labels=hc$labels[hc$order]
cut=cutree(hc,8)
#dend=color_branches(as.dendrogram(hc),k=length(unique(cut)),
#                    col=labelcolor[unique(cut[labels])])
dend <- as.dendrogram(hc)

circos.clear()
pdf(paste0(out_basepath, "plots/circos_barplot_nmr.perc.pdf"),width = 20, height = 20)
circos.par(cell.padding=c(0,0,0,0), start.degree = 45)
circos.initialize("a",xlim=c(0,nrow(t)))

circos.track(ylim=c(0,1),track.height=.2,track.margin=c(0,0),bg.border=NA,
             panel.fun=function(x,y)for(i in 1:nrow(t))circos.text(i-.5,0,labels[i],adj=c(0,.5),
                                                                   facing="clockwise",niceFacing=T,cex=.75,col= "black"))

circos.track(ylim=c(0,1),track.height=.3,track.margin=c(0,.01),bg.border=NA,
             panel.fun=function(x,y)circos.barplot(as.matrix(t)[hc$order,],-.5+1:nrow(t),
                                                   col=barcolor,bar_width=1,lwd=.3,border="gray20"))

circos.track(ylim=c(0,attr(dend,"height")),track.height=.4,track.margin=c(0,0),
             bg.border=NA,panel.fun=function(x,y)circos.dendrogram(dend))

circos.clear()
dev.off()

# 
# d <- read.delim("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/results/data/LLD_bloodlipids_nmr.txt", header = T, sep = "\t", as.is = T, check.names = F, row.names = 1)
# d <- na.omit(d)
# cormat <- cor(d, method = "spearman")
# km <- kmeans(cor2dist(cormat), num_k, nstart = 50)
# cl <- km$cluster
# clus <- cl[match(row.names(res_rsq_table), names(cl))]
# ord <- order(clus, res_rsq_table$sex)
# 
# # grouped <- res_rsq_table
# # grouped$group <- 0
# # grouped[grouped$age > 0.1 && grouped$sex > 0.1, "group"] <- 3
# # grouped[grouped$age < 0.1 && grouped$sex > 0.1, "group"] <- 2
# # grouped[grouped$age > 0.1 && grouped$sex < 0.1, "group"] <- 1
# # ord <-row.names(grouped[order( grouped$group),])
# 
# res_rsq_table_long <- melt(as.matrix(res_rsq_table[ord,]))
# colnames(res_rsq_table_long) <- c("phenotype", "model", "explained_variance")
# res_rsq_table_long$model <- factor(res_rsq_table_long$model, levels = c("inter", "age", "sex","null_model"))
# res_rsq_table_long$phenotype <- factor(res_rsq_table_long$phenotype, levels = unique(res_rsq_table_long$phenotype))
# 
# pdf(paste0(out_basepath, "plots/NMR_exp_variance3.pdf"), width = 5, height = 30)
# 
# ggplot(res_rsq_table_long, aes(width = 1, fill = model, y = explained_variance, x = phenotype)) + 
#   geom_bar(position="stack", stat="identity") +
#   scale_fill_brewer(palette = "Set2") +
#   theme_minimal() + coord_flip()
# dev.off()
