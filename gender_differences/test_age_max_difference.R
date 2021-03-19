library(RColorBrewer)

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

d <- read.delim("../all_LL/ttest_res.txt", sep = "\t", header = T, as.is = T, check.names = F)
d$w_p_log <- -log10(d$w_p)
d$m_p_log <- -log10(d$m_p)

pdf("../all_LL/ttest_max_differences2.pdf", height = 5, width = 10)
par(mfrow=c(1,2))
palette(c(col2transparent("indianred1", 130),col2transparent("dodgerblue1", 130)))
plot(m_p_log~age, data = d, pch = 16, frame.plot = F, axes = T, col = 2, ylim = c(0, max(d$w_p_log, d$m_p_log)), ylab = "-log10(P)", main = "men")
abline(h = pretty(d$w_p_log), col = "grey")
abline(v = pretty(d$age), col = "grey")

plot(w_p_log~age, data = d, pch = 16, frame.plot = F, axes = T, col = 1, ylim = c(0, max(d$w_p_log, d$m_p_log)), ylab = "-log10(P)", main = "women")
abline(h = pretty(d$w_p_log), col = "grey")
abline(v = pretty(d$age), col = "grey")
dev.off()

pdf("../all_LL/ttest_max_differences_coloured.pdf", height = 5, width = 10)
par(mfrow=c(1,2))
palette(brewer.pal(n = 12, name = "Set3"))
plot(m_p_log~age, data = d, pch = 16, frame.plot = F, axes = T, col = as.factor(pheno), ylim = c(0, max(d$w_p_log, d$m_p_log)), ylab = "-log10(P)", main = "men")
abline(h = pretty(d$w_p_log), col = "grey")
abline(v = pretty(d$age), col = "grey")

plot(w_p_log~age, data = d, pch = 16, frame.plot = F, axes = T, col = as.factor(pheno), ylim = c(0, max(d$w_p_log, d$m_p_log)), ylab = "-log10(P)", main = "women")
abline(h = pretty(d$w_p_log), col = "grey")
abline(v = pretty(d$age), col = "grey")
dev.off()


pmax_w <- with(d, ave(w_p_log, pheno, FUN=function(x) seq_along(x)==which.max(x)))==1
d_pmax_w <- d[pmax_w, ]
pmax_m <- with(d, ave(m_p_log, pheno, FUN=function(x) seq_along(x)==which.max(x)))==1
d_pmax_m <- d[pmax_m, ]
pdf("../all_LL/ttest_max_differences_minP.pdf", height = 5, width = 10)
par(mfrow=c(1,2))
palette(c(col2transparent("indianred1", 130),col2transparent("dodgerblue1", 130)))
plot(m_p_log~age, data = d_pmax_m, pch = 16, frame.plot = F, axes = T, col = 2, ylim = c(0, max(d$w_p_log, d$m_p_log)), ylab = "-log10(P)", main = "men")
abline(h = pretty(d$w_p_log), col = "grey")
abline(v = pretty(d$age), col = "grey")

plot(w_p_log~age, data = d_pmax_w, pch = 16, frame.plot = F, axes = T, col = 1, ylim = c(0, max(d$w_p_log, d$m_p_log)), ylab = "-log10(P)", main = "women")
abline(h = pretty(d$w_p_log), col = "grey")
abline(v = pretty(d$age), col = "grey")
dev.off()