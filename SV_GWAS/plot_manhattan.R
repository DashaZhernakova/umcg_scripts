library(ggplot2)
library(tidyr)
library(dplyr)
#library(ggrepel)

setwd("/data/umcg-tifn/SV/SV_GWAS/results_all_summary_stats/dSV/meta_combined_fdr")
fname = "eQTLs.txt.gz"
d <- read.table(gzfile(fname), header = T, as.is = T, sep = "\t", check.names = F)

d$CHR <- as.numeric(d[,3])
d$BP <- as.numeric(d[,4])

d$bac <- d[,5]
d$P <- d[,1]


don <- d %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(d, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)
  
  # Add highlight and annotation information
  #mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  #mutate( is_annotate=ifelse(-log10(P)>8, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
#pdf("dSV_manhattan_noannot.v5.pdf",  height =8, width = 20, useDingbats = F)
png("/data/umcg-tifn/SV/SV_GWAS/plots/dSV_manhattan_noannot.v5.png", type = "cairo", width = 15, height = 5, units = 'in', res = 400)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
#ggplot(don, aes(x=BPcum, y=log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "#9ECAE1"), 22 )) +
  geom_hline(yintercept=8.302771, color = "#EF3B2C") +
  #geom_hline(yintercept=-9.484126, color = "#EF3B2C", size = 0.3) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) , limits = c(0,45)) +     # remove space between plot area and x axis
  #scale_y_continuous(expand = c(0, 0) , limits = c(-45,0)) +     # remove space between plot area and x axis
  
  # Add highlighted points
  #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=bac), size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()
