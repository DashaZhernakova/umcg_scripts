library(ggplot2)

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

d <- read.delim("../all_LL/ttest_res.txt", sep = "\t", header = T, as.is = T, check.names = F)
d$w_p_log <- -log10(d$w_p)
d$m_p_log <- -log10(d$m_p)

pdf("../ttest_max_differences.pdf", height = 5, width = 5)
palette(c(col2transparent("indianred1", 130),col2transparent("dodgerblue1", 130)))
plot(m_p_log~age, data = d, pch = 16, frame.plot = F, axes = T, col = 2, ylim = c(0, max(d$w_p_log, d$m_p_log)))
abline(h = pretty(d$w_p_log), col = "grey")
abline(v = pretty(d$age), col = "grey")

points(w_p_log~age, data = d, pch = 16, col = 1)
dev.off()