library('mgcv')
library('dplyr')
setwd("C:/Users/Dasha/work/UMCG/data/gender_differences/omics/plots_all_pheno/diseases/")
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/preprocessing_gam_fitting_functions.R")
source("C:/Users/Dasha/work/UMCG/umcg_scripts/gender_differences/additional_functions.R")
pheno_path <- "../../data/LLD_covariates_291220.txt"
traits_path <- "C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/LLD_diseases_1660.txt"


col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

merge_match <- function(d1,d2){
  d1_m <- d1[match(row.names(d2), row.names(d1), nomatch = 0 ),]
  d2_m <- d2[match(row.names(d1_m), row.names(d2), nomatch = 0),]
  print(all(row.names(d1_m) == row.names(d2_m)))
  merged <- cbind(d1_m, d2_m)
  colnames(merged) <- c(colnames(d1_m), colnames(d2_m))
  return (merged)
}

# traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- as.data.frame(sapply(traits0, function(x) as.numeric(as.character(x))))
row.names(traits) <- row.names(traits0)

# Age, gender and other phenotypes
pheno0 <- read.delim(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- na.omit(pheno0)

traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))
m <- cbind(pheno_m[,c("age","gender_F1M2")], traits_m)
min_age <- 20
max_age <- 80
m <- m[m$age >= min_age & m$age < max_age,]

prev_w <- aggregate(m[m$gender_F1M2 == 1,], by = list(m[m$gender_F1M2 == 1,"age"]), FUN = mean, na.rm = TRUE)
prev_m <- aggregate(m[m$gender_F1M2 == 2,], by = list(m[m$gender_F1M2 == 2,"age"]), FUN = mean,  na.rm = TRUE)
rownames(prev_w) <- prev_w$age
rownames(prev_m) <- prev_m$age

pheno_idx = 40

###

tmp_w <- prev_w[, c(1, pheno_idx)]
tmp_m <- prev_m[, c(1, pheno_idx)]
colnames(tmp_w) <- c("age", "phen")
colnames(tmp_m) <- c("age", "phen")

rownames(tmp_w) <- tmp_w$age
rownames(tmp_m) <- tmp_m$age
x <- merge_match(tmp_w, tmp_m)
x <- x[c(2,4)]
colnames(x) <- c("wom", "men")
x2 <- t(x)
bp <- barplot(as.matrix(x2), beside = T, border="white", col=c("indianred1","dodgerblue1"),
              xlab = "age", ylab = "disease prevalence", main = paste0(colnames(prev_m)[pheno_idx]))




###
idx = pheno_idx - 3
trait_name = colnames(traits_m)[idx]
outlier_correction_method = "zscore"
log_transform = F
scale_transform = F
merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, n_points = 300, make_plots = T, gam_family = binomial(link = "logit"), label = '', add_breakpoints = F)
merged_tab <- merged_tab[merged_tab$age > min_age & merged_tab$age <= max_age,]

##

hist(merged_tab[merged_tab$gender_F1M2 == 1 & merged_tab[,1] == 1, "age"], breaks = c(seq(min_age,max_age,5)),col = col2transparent("indianred1",70) , xlim = c(min_age,max_age), main = paste0(trait_name , " = 1"), xlab = "age")
hist(merged_tab[merged_tab$gender_F1M2 == 2 & merged_tab[,1] == 1, "age"], breaks = c(seq(min_age,max_age,5)), col = col2transparent("dodgerblue1", 70), xlim = c(min_age,max_age), add = T)

hist(merged_tab[merged_tab$gender_F1M2 == 2 & merged_tab[,1] == 0, "age"], breaks = c(seq(min_age,max_age,5)), col = col2transparent("dodgerblue1", 70), xlim = c(min_age,max_age), main = paste0(trait_name , " = 0"), xlab = "age")
hist(merged_tab[merged_tab$gender_F1M2 == 1 & merged_tab[,1] == 0, "age"], breaks = c(seq(min_age,max_age,5)),col = col2transparent("indianred1",70) , xlim = c(min_age,max_age), add = T)

draw_hist(merged_tab)

#
### Adapted from https://rdrr.io/cran/popbio/src/R/logi.hist.plot.R
#
draw_hist <- function(merged_tab, scale.hist = 5, las.h1 = 1, ylabel2 = "Frequency"){
  
  col.hist.w = col2transparent("indianred1",70)
  col.hist.m = col2transparent("dodgerblue1", 70)
  h.br.step <- 5
  h.br <- seq(min_age, max_age, h.br.step)
  
  independ.w<- merged_tab[merged_tab$gender_F1M2 == 1, "age"]
  depend.w <- merged_tab[merged_tab$gender_F1M2 == 1, 1]
  independ.m <- merged_tab[merged_tab$gender_F1M2 == 2, "age"]
  depend.m <- merged_tab[merged_tab$gender_F1M2 == 2, 1]
  
  
  h.x <- hist(independ[depend == 0],
              breaks = h.br,
              plot = FALSE
  )$mid
  #
  # get counts in each bin
  #
  #men
  h.m.y0 <- hist(independ.m[depend.m == 0],
               breaks = h.br,
               plot = FALSE
  )$counts
  h.m.y1 <- hist(independ.m[depend.m == 1],
               breaks = h.br,
               plot = FALSE
  )$counts
  # women
  h.w.y0 <- hist(independ.w[depend.w == 0],
                 breaks = h.br,
                 plot = FALSE
  )$counts
  h.w.y1 <- hist(independ.w[depend.w == 1],
                 breaks = h.br,
                 plot = FALSE
  )$counts
  
  
  # scale the histogram bars to max desired length:
  scale_val <- max(c(h.m.y0, h.m.y1, h.w.y0, h.w.y1)) * scale.hist
  
  h.m.y0n <- h.m.y0 / scale_val
  h.m.y1n <- 1 - h.m.y1 / scale_val
  h.w.y0n <- h.w.y0 / scale_val
  h.w.y1n <- 1 - h.w.y1 / scale_val
  
  # draw bottom histogram:
  for (i in 1:length(h.m.y0n)) {
    if (h.m.y0n[i] > 0) {
      polygon(c(rep(h.br[i], 2), rep(h.br[i] + h.br.step/2, 2)),
              c(0, rep(h.m.y0n[i], 2), 0),
              col = col.hist.m
      )
    }
    if (h.w.y0n[i] > 0) {
      polygon(c(rep(h.br[i] + h.br.step/2, 2), rep(h.br[i + 1], 2)),
              c(0, rep(h.w.y0n[i], 2), 0),
              col = col.hist.w
      )
    }
  }
  # draw top histogram:
  for (i in 1:length(h.y1n)) {
    if (h.m.y1n[i] < 1) {
      polygon(c(rep(h.br[i], 2), rep(h.br[i] + h.br.step/2, 2)),
              c(h.m.y1n[i], 1, 1, h.m.y1n[i]),
              col = col.hist.m
      )
    }
    if (h.w.y1n[i] < 1) {
      polygon(c(rep(h.br[i] + h.br.step/2, 2), rep(h.br[i + 1], 2)),
              c(h.w.y1n[i], 1, 1, h.w.y1n[i]),
              col = col.hist.w
      )
    }
  }
  
  axis.hist <- function(h.m.y0, h.m.y1, h.w.y0, h.w.y1, scale.hist,
                        las = las.h1) {
    tope <- max(c(h.m.y0, h.m.y1, h.w.y0, h.w.y1))
    label.down <- c(
      0, (ceiling(tope / 10)) * 5,
      (ceiling(tope / 10)) * 10
    )
    label.up <- c(
      (ceiling(tope / 10)) * 10,
      (ceiling(tope / 10)) * 5, 0
    )
    at.down <- label.down / (tope * scale.hist)
    at.up <- 1 - (label.up / (tope * scale.hist))
    at.hist <- c(at.down, at.up)
    label.hist <- c(label.down, label.up)
    axis(
      side = 4, at = at.hist, labels = label.hist,
      las = las
    )
    #mtext(ylabel2, side = 4, line = 3, cex = 0.6)
    text(par("usr")[2] * 1.05,0.6, ylabel2, srt = -90, xpd = TRUE, pos = 4)
  }
  axis.hist(h.m.y0, h.m.y1, h.w.y0, h.w.y1, scale.hist)
  axis(side = 2, las = las.h1)
}

draw_disease_prevalence <- function(merged_tab){
  prevalence_w <- aggregate(merged_tab[merged_tab$gender_F1M2 == 1, 1], by = list(merged_tab[merged_tab$gender_F1M2 == 1,"age"]), FUN = mean, na.rm = TRUE)
  prevalence_m <- aggregate(merged_tab[merged_tab$gender_F1M2 == 2, 1], by = list(merged_tab[merged_tab$gender_F1M2 == 2,"age"]), FUN = mean,  na.rm = TRUE)
  
  colnames(prevalence_w) <- c("age", "phen")
  colnames(prevalence_m) <- c("age", "phen")
  
  rownames(prevalence_w) <- prevalence_w$age
  rownames(prevalence_m) <- prevalence_m$age
  x <- merge_match(prevalence_w, prevalence_m)
  x <- x[c(2,4)]
  colnames(x) <- c("wom", "men")
  x2 <- t(x)
  bp <- barplot(as.matrix(x2), beside = T, border="white", col=c("indianred1","dodgerblue1"),
                xlab = "age", ylab = "disease prevalence", main = paste0(colnames(prev_m)[pheno_idx]))
  
}