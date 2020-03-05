library(reshape2)
library(ggplot2)
library(pheatmap)
res_dif_all <- read.delim("../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.txt", header = T, sep = "\t", as.is = T, check.names = F)
res_dif_all <- subset(res_dif_all, select = -age)
cor_matrix <- cor(res_dif_all)
cor_matrix <- reorder_cormat(cor_matrix)

pdf('../anova+gam/rnaseq_plots/ageDEgenes/LLD_expression_ageDEgenes_result.diff.heatmap.pdf', width = 15, height = 15, useDingbats = F)
pheatmap(cor_matrix)
dev.off()




# ggplot version
cor_matrix_mlt <- melt(cor_matrix)

ggheatmap <- ggplot(cor_matrix_mlt, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

reorder_cormat <- function(cormat){
 # Use correlation between variables as distance
 dd <- as.dist((1-cormat)/2)
 hc <- hclust(dd)
 cormat <-cormat[hc$order, hc$order]
}