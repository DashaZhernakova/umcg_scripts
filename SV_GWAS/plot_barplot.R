d <- read.delim("dSV_ABO_A_counts.txt", sep = "\t", as.is = T, check.names = F)
library(ggplot2)
d$abo_fut2 <- interaction(d$FUT2_status, d$ABO_A)
d$abo_fut2 <- factor(d$abo_fut2, levels = c("secretor.A", "secretor.O/B", "non-secretor.A", "non-secretor.O/B"))
d$cohort <- as.factor(d$cohort)
d$cohort <- factor(d$cohort, levels = c("DAG3", "LLD", "500FG"))

pdf("abo_a_vs_dSV.pdf")
ggplot(d[d$sv == 2,], aes (fill = cohort, y = frac, x = abo_fut2)) +
  geom_bar(position="dodge", stat="identity") + 
  theme_classic() +
  scale_fill_viridis_d(alpha = 0.5) +
  theme(text = element_text(size = 12), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.6)) +
  labs(y="Fraction of samples without deletion")
dev.off()
