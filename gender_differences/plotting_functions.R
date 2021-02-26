library(RColorBrewer)
library('dplyr')

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

draw_plot <- function(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot = T, breakpoints = NULL, factor_name = "", alpha_points = 40, breakpoints_intervals = NULL, label = ""){
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
  binaryPhenotype <- F
  if (all(ylims == c(0,1))) {
    ylabel <- "Probability"
    binaryPhenotype <- T
  }
  
  
  ## draw base plot
  #palette(c(col2transparent("indianred1", alpha_points),col2transparent("dodgerblue1", alpha_points)))
  palette(c(col2transparent("#ff9999", 120),col2transparent("#99ccff", 120)))
  par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i") # Remove plot padding
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  if (add_inter_p_to_plot) {
    plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction p = ", format(gam.p, digits = 3))
  } else {
    plot_title <- pheno_name
  }
  if (! binaryPhenotype){
    plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
         main = plot_title, 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
         ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
         xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
    
    if (length(factor_name) > 0){
      points(merged_tab2[merged_tab2[,factor_name] == 1,"age"], merged_tab2[merged_tab2[,factor_name],"phenotype"], pch = 8, col = col2transparent("gold", 120), cex = 0.6)
    }
  } else { # do not add points for binary pheno
    plot(1, type="n",
         main = paste0(pheno_name, '\n', label, "\nGAM interaction p = ", format(gam.p, digits = 3)), 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
         ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
         xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
    
    add_prevalence_hist_to_plot(merged_tab2)
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
