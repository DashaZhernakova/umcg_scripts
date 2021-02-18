library(RColorBrewer)
library('dplyr')

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

draw_plot <- function(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, breakpoints = F, factor_name = "", alpha_points = 40, breakpoints_intervals = NULL, label = ""){
  cex_main = 1
  ylims <- with(merged_tab, range(phenotype))
  ylabel <- pheno_name
  if (nchar(pheno_name) > 40){
    spl <- unlist(strsplit(pheno_name, "\\|"))
    ylabel <- spl[length(spl)]
    cex_main <- 0.8
    if (nchar(pheno_name) > 50) cex_main <- 0.7
    if (nchar(pheno_name) > 60) cex_main <- 0.6
    
  }
  if (all(ylims == c(0,1))) ylabel <- "Probability"
  
  
  ## draw base plot
  #palette(c(col2transparent("indianred1", alpha_points),col2transparent("dodgerblue1", alpha_points)))
  palette(c(col2transparent("#ff9999", 120),col2transparent("#99ccff", 120)))
  par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i") # Remove plot padding
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
       main = paste0(pheno_name, '\n', label, "\nGAM interaction p = ", format(gam.p, digits = 3)), 
       cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
       ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
       xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
  
  if (length(factor_name) > 0){
    points(merged_tab2[merged_tab2[,factor_name] == 1,"age"], merged_tab2[merged_tab2[,factor_name],"phenotype"], pch = 8, col = col2transparent("gold", 120), cex = 0.6)
  }
  
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2$age), col = "grey90")
  
  #at = pretty(merged_tab2$age)
  #mtext(side = 1, text = at, at = at, 
  #      col = "grey20", line = 1, cex = 0.4)
  
  #at = pretty(merged_tab2$phenotype)  
  #mtext(side = 2, text = at, at = at, col = "grey20", line = 1, cex = 0.4)
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
    polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(cols[[l]], 65), border = NA)
  }
  
  if (all(ylims == c(0,1))) add_prevalence_hist_to_plot(merged_tab2)
  
  if(length(breakpoints) > 0){
    br_w <- breakpoints[[1]]
    br_m <- breakpoints[[2]]
    for (br in br_w){
      abline(v = as.numeric(br), col = cols[1], lty = 2)
    }
    for (br in br_m){
      abline(v = as.numeric(br), col = cols[2], lty = 2)
    }
  }
  ymin <- min(pretty(merged_tab2$phenotype), min(merged_tab2$phenotype) - 1)
  ymax <- max(pretty(merged_tab2$phenotype), max(merged_tab2$phenotype) + 2)
  if(!is.null(breakpoints_intervals)){
    br_w <- breakpoints_intervals[[1]]
    br_m <- breakpoints_intervals[[2]]
    for (br in br_m){
      rect(br[1], ymin, br[2], ymax, col = 2, border = 2, density = -1)
    }
    for (br in br_w){
      rect(br[1], ymin, br[2], ymax, col = 1, border = 1, density = -1)
    }
  }
}

# plot more than two fitted lines 
draw_plot_multiline <- function(merged_tab, pheno_name, pdat_list, color_list, factor_name = "", min_age = 20, max_age = 80, label = ''){
  cex_main = 1
  ylims <- with(merged_tab, range(phenotype))
  ylabel <- pheno_name
  if (nchar(pheno_name) > 40){
    spl <- unlist(strsplit(pheno_name, "\\|"))
    ylabel <- spl[length(spl)]
    cex_main <- 0.8
    if (nchar(pheno_name) > 50) cex_main <- 0.7
    if (nchar(pheno_name) > 60) cex_main <- 0.6
    
  }
  if (all(ylims == c(0,1))) ylabel <- "Probability"
  
  
  ## draw base plot
  palette(c(col2transparent("#ff9999", 120),col2transparent("#99ccff", 120)))
  par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i") # Remove plot padding
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
       main = pheno_name, 
       cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
       ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
       xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
  
  if (factor_name != ""){
    subs <- merged_tab2[merged_tab2[,factor_name] == 1,]
    points(subs$age, subs$phenotype, pch = 8, col = col2transparent("gold", 120), cex = 0.6)
  }
  
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2$age), col = "grey90")
  
  ## add the fitted lines
  for (i in 1:length(pdat_list)) {
    dd <- pdat_list[[i]]
    lines(pred ~ age, data = dd, col = color_list[[i]], lwd = 2)
    polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(color_list[[i]], 65), border = NA)
  }
}

calculate_sex_diff_ttest <- function(merged_tab, covariates = c(), min_age = 20, max_age = 80){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  
  if (length(covariates) == 0){
    resid <- merged_tab$phenotype 
    #res <- lm(phenotype ~ gender_F1M2, data = merged_tab)
  } else{
    resid <- residuals( lm(as.formula(paste("phenotype ~ ",  paste(covariates, collapse = "+"))), data = merged_tab))
  }
  #sex_p <- summary(res)$coefficients[2,4]
  t_res <- t.test(resid[merged_tab$gender_F1M2 == 1], resid[merged_tab$gender_F1M2 == 2])
  sex_p <- t_res$p.value
  return(sex_p)
}

draw_plot0 <- function(merged_tab, pheno_name, pdat){
  cex_main = 1
  ylims <- with(merged_tab, range(phenotype))
  ylabel <- pheno_name
  if (nchar(pheno_name) > 40){
    spl <- unlist(strsplit(pheno_name, "\\|"))
    ylabel <- spl[length(spl)]
    cex_main <- 0.8
    if (nchar(pheno_name) > 50) cex_main <- 0.7
    if (nchar(pheno_name) > 60) cex_main <- 0.6
    
  }
  
  
  ## draw base plot
  palette(c(col2transparent("indianred1", 40),col2transparent("dodgerblue1", 40)))
  plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
       main = paste0(pheno_name, ' ', label, "\nGAM interaction p = ", format(gam.p, digits = 3)), 
       cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main)
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
    polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(cols[[l]], 50), border = NA)
  }
  
  breakpoints <- NULL
  if (add_breakpoints){
    breakpoints <- get_breakpoints(merged_tab)
    br_w <- breakpoints[[1]]
    br_m <- breakpoints[[2]]
    for (br in br_w){
      abline(v = as.numeric(br), col = cols[1], lty = 2)
    }
    for (br in br_m){
      abline(v = as.numeric(br), col = cols[2], lty = 2)
    }
  }
  
  #if (! is.null(res_dif = NULL)){
  #  #Plot the difference
  #  plot(res_dif$age, res_dif$diff, type = 'l')
  #}
}

test_polynomial_interaction <- function(merged_tab, model_fam = NO(), gm_mean_cov = F){
  if (! gm_mean_cov){
    #lgam.fit <- gamlss(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) + I(age^4), data = merged_tab, family = model_fam, method = RS(30))
    #lm2 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^3) + I(age^4) + gender_F1M2:(age + I(age^2) + I(age^3) + I(age^4)), data = merged_tab, family = model_fam, method = RS(30)) 
    lgam.fit <- gamlss(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) , data = merged_tab, family = model_fam, method = RS(30))
    lm2 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)), data = merged_tab, family = model_fam, method = RS(30)) 
  } else {
    lgam.fit <- gamlss(phenotype ~ gender_F1M2 + age + gm_mean + I(age^2) + I(age^3) , data = merged_tab, family = model_fam, method = RS(30))
    lm2 <- gamlss(phenotype ~ gender_F1M2 + age + gm_mean + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)), data = merged_tab, family = model_fam, method = RS(30)) 
    
  }
  
  lr <- LR.test(lgam.fit, lm2, print = F)
  m <- lgam.fit
  if (lr$p.val < 0.05){
    m <- lm2
  }
  return (list("lr_p.val" = lr$p.val, "mod" = m))
}


# Calculates Cohen's R^2 to estimate the local effect size of the variable of interest
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3328081/
# According to Cohen’s (1988) guidelines, f2>=0.02, f2>=0.15, and f2>=0.35 represent small, medium, and large effect sizes, respectively.
calculate_cohens_f2 <- function(mod_full, mod_part){
  r_full <- summary(mod_full)$r.sq
  r_sub <- summary(mod_part)$r.sq
  f2 <- (r_full - r_sub)/(1-r_full)
  return(f2)
}

correct_for_covariates_before <- function(traits_m, pheno_m, covariates){
  new_traits <- matrix(nrow = nrow(traits_m), ncol = ncol(traits_m))
  colnames(new_traits) <- colnames(traits_m)
  rownames(new_traits) <- rownames(traits_m)
  
  if (! all(covariates %in% colnames(pheno_m))){
    print("Error! Some covariates are not present in the phenotype file")
    print(paste0("Covariates: ", paste(covariates, collapse = ",")))
    print(paste0("Covariate file phenotypes: ", paste(colnames(pheno_m), collapse = ",")))
  }
  
  for (col in 1:ncol(traits_m)){
    d <- cbind(traits_m[,col], pheno_m[,covariates])
    colnames(d)[1] <- "phenotype"
    lm.fit <- lm(as.formula(paste("phenotype ~  ", paste(covariates, collapse = "+"), sep = " ")), data = d)
    new_traits[,col] <- residuals(lm.fit)
  }
  return(new_traits)
}

#
### Adapted from https://rdrr.io/cran/popbio/src/R/logi.hist.plot.R
#
add_hist_to_plot <- function(merged_tab, scale.hist = 5, las.h1 = 1, ylabel2 = "Count"){
  
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
    text(par("usr")[2] * 1.08,0.6, ylabel2, srt = -90, xpd = TRUE, pos = 4)
  }
  axis.hist(h.m.y0, h.m.y1, h.w.y0, h.w.y1, scale.hist)
  axis(side = 2, las = las.h1)
}

#
### Adapted from https://rdrr.io/cran/popbio/src/R/logi.hist.plot.R
#
add_prevalence_hist_to_plot <- function(merged_tab, scale.hist = 5, las.h1 = 1, ylabel2 = "Count"){
  
  col.hist.w = col2transparent("indianred1",70)
  col.hist.m = col2transparent("dodgerblue1", 70)
  h.br.step <- 3
  h.br <- seq(min_age, max_age, h.br.step)
  
  independ.w<- merged_tab[merged_tab$gender_F1M2 == 1, "age"]
  depend.w <- merged_tab[merged_tab$gender_F1M2 == 1, 1]
  independ.m <- merged_tab[merged_tab$gender_F1M2 == 2, "age"]
  depend.m <- merged_tab[merged_tab$gender_F1M2 == 2, 1]
  
  
  h.x <- hist(independ.m[depend.m == 0],
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
  
  
   
  h.m <- h.m.y1 / (h.m.y1 + h.m.y0)
  h.w <- h.w.y1 / (h.w.y1 + h.w.y0)
  
  # draw histogram:
  for (i in 1:length(h.m)) {
    if (h.m[i] > 0) {
      polygon(c(rep(h.br[i], 2), rep(h.br[i] + h.br.step/2, 2)),
              c(0, rep(h.m[i], 2), 0),
              col = col.hist.m, border = "darkgrey"
      )
    }
    if (h.w[i] > 0) {
      polygon(c(rep(h.br[i] + h.br.step/2, 2), rep(h.br[i + 1], 2)),
              c(0, rep(h.w[i], 2), 0),
              col = col.hist.w, border = "darkgrey"
      )
    }
  }
  
  axis(
    side = 4, at = seq(0,1,0.2), labels = seq(0,1,0.2),
    las = 1
  )
  text(par("usr")[2] * 1.08, 0.75, "Frequency", srt = -90, xpd = TRUE, pos = 4)
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

run_for_split_by_covariate <- function(merged_tab, pheno_name, covariate_to_split, highlight_positive = F, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots = T,  gam_family = gaussian(), min_age = 20, max_age = 80){
  pdat_list <- list()
  col_list <- list()
  colors <- list("indianred1", "orange1", "dodgerblue1", "cadetblue1")
  w0 <- merged_tab[merged_tab$gender_F1M2 == 1 & merged_tab[,covariate_to_split] == 0,]
  w1 <- merged_tab[merged_tab$gender_F1M2 == 1 & merged_tab[,covariate_to_split] == 1,]
  m0 <- merged_tab[merged_tab$gender_F1M2 == 2 & merged_tab[,covariate_to_split] == 0,]
  m1 <- merged_tab[merged_tab$gender_F1M2 == 2 & merged_tab[,covariate_to_split] == 1,]  
  
  cnt <- 1
  cnt2 <- 1
  for (dat in list(w0,w1,m0,m1)){
    if (nrow(dat) > 0){
      pdat_list[[cnt]] <- fit_gam_age_only(dat, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots = T,  gam_family = gaussian(), min_age = 20, max_age = 80)
      col_list[[cnt]] <- colors[[cnt2]]
      cnt <- cnt + 1
    }
    cnt2 <- cnt2 + 1
  }
  factor_to_highlight <- ""
  if (highlight_positive) factor_to_highlight <- covariate_to_split
  draw_plot_multiline(merged_tab, pheno_name, pdat_list, col_list, factor_name = factor_to_highlight)
}
