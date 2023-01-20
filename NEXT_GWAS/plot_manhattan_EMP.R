library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)

folder = args[1]
d <- read.delim(gzfile(paste0("/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/", folder, "/eQTLs.txt.gz")), header = T, as.is = T, sep = "\t", check.names = F)
signif <- read.delim(paste0("/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/", folder, "/eQTLsFDR0.05-ProbeLevel.txt"), header = T, as.is = T, sep = "\t", check.names = F)
#top_res <- read.delim(paste0("/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/", folder, "/eQTLProbesFDR0.05-ProbeLevel.txt"), header = T, as.is = T, sep = "\t", check.names = F)

#setwd("/Users/Dasha/work/UMCG/data/NEXT/HMO_GWAS/EMP_HMO_M2/")
#d <- read.delim(gzfile("eQTLs.txt.gz"), header = T, as.is = T, sep = "\t", check.names = F)
#signif <- read.delim("eQTLsFDR0.05-ProbeLevel.txt", header = T, as.is = T, sep = "\t", check.names = F)
#top_res <- read.delim("eQTLProbesFDR0.05-ProbeLevel.txt", header = T, as.is = T, sep = "\t", check.names = F)

fdr_cutoff <- max(signif$PValue)
d$SNP_pheno <- paste0(d$SNPName, "-", d$ProbeName)
signif$SNP_pheno <- paste0(signif$SNPName, "-", signif$ProbeName)
to_annotate <- c()

for (chr in seq(1,22)){
  for (pheno in unique(signif$ProbeName)){
    subs <- signif[signif$ProbeName == pheno & signif$SNPChr == as.character(chr),]
    if (nrow(subs) > 0) to_annotate <- c(to_annotate, subs[1,"SNP_pheno"])
    #print(to_annotate)
  }
}


d$is_annotate <- ifelse(d$SNP_pheno %in% to_annotate, "yes", "no")

d$SNPChrPos <- as.numeric(d$SNPChrPos)


don <- d %>% 
  
  # Compute chromosome size
  group_by(SNPChr) %>% 
  summarise(chr_len=max(SNPChrPos)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(d, ., by=c("SNPChr"="SNPChr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(SNPChr, SNPChrPos) %>%
  mutate( BPcum=SNPChrPos+tot)
  
  # Add highlight 
  #mutate( is_highlight=ifelse(SNP %in% snps_to_highlight, "yes", "no") %>%


axisdf <- don %>% group_by(SNPChr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )


png(paste0("/groups/umcg-llnext/tmp01/umcg-dzhernakova/HMO_GWAS/plots/", folder, "_manhattan_noannot.png"), width = 15, height = 5, units = 'in', res = 400)
ggplot(don, aes(x=BPcum, y=-log10(PValue))) +
  # Show all points
  geom_point( aes(color=as.factor(SNPChr)),  size=1) +
  scale_color_manual(values = rep(c("#495DA0", "#74AFDF"), 22 )) +
  geom_hline(yintercept=8.30103, color = "#EF3B2C") +
  geom_hline(yintercept=-1*log(fdr_cutoff,10), color = "orange") +
  # custom X axis:
  xlab("Chromosome") +
  scale_x_continuous( label = axisdf$SNPChr, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) , limits = c(0,45)) +
  
  # Add annotation
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=ProbeName), size=2, max.overlaps = 20) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line.y = element_line(color="lightgrey", size = 0.5),
    axis.ticks.x = element_blank()
    
  )
dev.off()
