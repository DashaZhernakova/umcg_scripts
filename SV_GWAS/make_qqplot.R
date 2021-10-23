library(ggplot2)
library(patchwork)
library(qqplotr)
library(dplyr)

inflation <- function(ps) {
  chisq <- qchisq(1 - ps, 1)
  lambda <- median(chisq) / qchisq(0.5, 1)
  lambda
}


setwd('/groups/umcg-lld/tmp01/other-users/umcg-dzhernakova/SV_GWAS')
svs=c("Faecalibacterium__cf.__prausnitzii__KLE1255:101", "Eubacterium__rectale__DSM__17629:176", 
      "Bifidobacterium__longum:19", "Blautia__obeum__ATCC__29174:117", 
      "Coprococcus__comes__ATCC__27758:13", "Blautia__wexlerae__DSM__19850:266")

plots = list()
cnt = 1
for (sv in svs){
    d <- read.delim(gzfile(paste0("results/meta/dSV_meta_", sv,"_1.tbl.gz")), header = T, sep = "\t", as.is = T, check.names = F)
    d <- rename(d, c("P" = "P-value"))
    qq <- ggplot(data = d, mapping = aes(sample = P)) +
      stat_qq_band() +
      stat_qq_line() +
      stat_qq_point() +
      labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
      
      annotate(
        geom = "text",
        x = -Inf,
        y = Inf,
        hjust = -0.15,
        vjust = 1 + 0.15 * 3,
        label = sprintf("?? = %.2f", inflation(d$P)),
        size = 8
      ) 
    plots[[cnt]] <- qq
    cnt = cnt + 1
}
png("qqplots_dSV.png")
wrap_plots(plots[[1]], plots[[2]], plots[[3]], plots[[4]], plots[[5]], plots[[6]])
dev.off()

