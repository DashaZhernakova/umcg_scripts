library(RColorBrewer)
library('dplyr')
library(ggplot2)

col2transparent <- function(col, transparency){
  colRgb <- col2rgb(col)
  dodgerblueTransparent <- rgb(colRgb[1,1], colRgb[2,1], colRgb[3,1], transparency, names = NULL, maxColorValue = 255)
}

colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}

draw_plot <- function(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot = T, plot_title = NULL, plot_points = T,  factor_name = character(0), alpha_points = 40,  ymax_hist = 1, label = "", ylims_usr = NULL){
  
  cex_main = 1
  ylims <- with(merged_tab, range(phenotype))

  #if (nchar(pheno_name) > 40){
    #spl <- unlist(strsplit(pheno_name, "\\|"))
    #ylabel <- spl[length(spl)]
    #cex_main <- 0.8
    #if (nchar(pheno_name) > 50) cex_main <- 0.7
    #if (nchar(pheno_name) > 60) cex_main <- 0.6
    
  #}
  pheno_name <- paste(strwrap(pheno_name, width = 40), collapse = "\n")
  ylabel <- pheno_name
  pheno_name <- gsub(" \\(.*", "", pheno_name)
  binaryPhenotype <- F
  if (all(ylims == c(0,1))) {
    ylabel <- paste0(pheno_name, " frequency")
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
  if (plot_title == ""){
    if (add_inter_p_to_plot) {
      if (gam.p == 0) {
	      plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P < 2.23e-308")
      } else {
	      plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P = ", format(gam.p, digits = 3))
      }
    } else {
      plot_title <- pheno_name
    }
  }
  
  if (! binaryPhenotype){
    if (plot_points){
      plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
           main = plot_title, 
           cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
           ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
           xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
      
      if (length(factor_name) > 0){
        points(merged_tab2[merged_tab2[,factor_name] == 1,"age"], merged_tab2[merged_tab2[,factor_name],"phenotype"], pch = 8, col = col2transparent("gold", 120), cex = 0.6)
      }
    } else{
      plot(1, type="n", 
           main = plot_title, 
           cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
           ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
           xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
    }
  } else { # do not add points for binary pheno
    cat(ymax_hist)
    plot(1, type="n",
         main = paste0(pheno_name, '\n', label, "\nGAM interaction p = ", format(gam.p, digits = 3)), 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
         ylim =c(0,ymax_hist),
         xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
    
    add_prevalence_hist_to_plot(merged_tab2)
  }
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2$age), col = "grey90")
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
    polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(cols[[l]], 100), border = NA)
  }
  
}

# plot more than two fitted lines 
draw_plot_multiline <- function(merged_tab, pheno_name, pdat_list, color_list, factor_name = "", min_age = 20, max_age = 80, n_points = 300, label = '',plot_points = T){
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
  if (all(ylims == c(0,1))) ylabel <- paste0(pheno_name, " frequency")
  pheno_name <- gsub(" \\(.*", "", pheno_name)
  
  ## draw base plot
  palette(c(col2transparent("#ff9999", 120),col2transparent("#99ccff", 120)))
  par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i") # Remove plot padding
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  if (plot_points){
    plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
         main = pheno_name, 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
         ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
         xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
    
    if (factor_name != ""){
      subs <- merged_tab2[merged_tab2[,factor_name] == 1,]
      points(subs$age, subs$phenotype, pch = 8, col = col2transparent("gold", 120), cex = 0.6)
    }
  } else{
    plot(1, type="n", 
         main = pheno_name, 
         cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
         ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
         xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
  }
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2$age), col = "grey90")
  
  ## add the fitted lines
  for (i in 1:length(pdat_list)) {
    dd <- pdat_list[[i]]
    lines(pred ~ age, data = dd, col = color_list[[i]], lwd = 2)
    polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(color_list[[i]], 50), border = NA)
  }
}

draw_multiple_fitted_lines <- function(fitted_matrix, sign_inters = colnames(fitted_matrix), min_age = 20, max_age = 80, n_points = 300, ylab = 'phenotype', plot_title = ''){
  
  # ## draw base plot
  # par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
  #     mgp = c(2, 0.4, 0), # Dist' plot to label
  #     las = 1, # Rotate y-axis text
  #     tck = -.01, # Reduce tick length
  #     xaxs = "i", yaxs = "i") # Remove plot padding
  
  plot(1, type="n",cex = 0.6, xlab = "age",  ylab = ylab, frame.plot = T, axes = T, main = plot_title, ylim = c(-1, 1),xlim = c(min_age,max_age))
#       ylim =c(min(fitted_matrix), max(fitted_matrix)),
       
  
  age_seq <- seq(min_age, max_age,length.out = n_points)
  
  ## add the fitted lines
  
  for (i in 1:ncol(fitted_matrix)) {
    
    ltype = 1
    lwd = 3
    cols <- c(col2transparent("indianred1", 120), col2transparent("dodgerblue1", 120))
    if (! colnames(fitted_matrix)[i] %in% sign_inters) {
      ltype = 2 
      lwd = 2
      cols <- c(col2transparent("indianred1", 90), col2transparent("dodgerblue1", 90))
    }
    lines(age_seq, fitted_matrix[1:n_points, i] , col = cols[1], lwd = lwd, lty = ltype)
    lines(age_seq, fitted_matrix[(n_points + 1):(n_points*2), i] , col = cols[2], lwd = lwd, lty = ltype)
    #text(x = 28, y = fitted_matrix[1, i], col = "indianred1", colnames(fitted_matrix)[i])
    #text(x = 28, y = fitted_matrix[301, i], col = "dodgerblue1", colnames(fitted_matrix)[i])
    #Sys.sleep(4)
  }
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
  
  #axis(
  #  side = 4, at = seq(0,1,0.2), labels = seq(0,1,0.2),
  #  las = 1
  #)
  #text(par("usr")[2] * 1.08, 0.75, "Frequency", srt = -90, xpd = TRUE, pos = 4)
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


draw_smooth_scatter <- function(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot = T, plot_title = NULL, label = ""){
  cex_main = 1
  cat("Plotting a density smooth plot")
  ylims <- with(merged_tab, range(phenotype))
  ylabel <- pheno_name
  if (nchar(pheno_name) > 40){
    spl <- unlist(strsplit(pheno_name, "\\|"))
    ylabel <- spl[length(spl)]
    cex_main <- 0.8
    if (nchar(pheno_name) > 50) cex_main <- 0.7
    if (nchar(pheno_name) > 60) cex_main <- 0.6
    
  }
  pheno_name <- gsub(" \\(.*", "", pheno_name)
  
  ## draw base plot
  par(mar = c(6, 6, 6, 3), # Dist' from plot to side of page
      mgp = c(2, 0.4, 0), # Dist' plot to label
      las = 1, # Rotate y-axis text
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i") # Remove plot padding
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  if (plot_title == ""){
    if (add_inter_p_to_plot) {
      if (gam.p == 0) {
        plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P < 2.23e-308")
      } else {
        plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P = ", format(gam.p, digits = 3))
      }
    } else {
      plot_title <- pheno_name
    }
  }
  
  #plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
  #     main = plot_title, 
  #     cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
  #     ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
  #     xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
  
  w_palette <- brewer.pal(n = 9, name = "Reds")
  m_palette <- brewer.pal(n = 9, name = "Blues")
  
  #
  # Women
  #
  w <- merged_tab[merged_tab2$gender_F1M2 == 1,]
  
  w_denscols <- densCols(w$age, w$phenotype, colramp=colorRampPalette(c("black", "white")))
  cols <-  colorRampAlpha(w_palette, n = 256, alpha = 0.4)
  w$dens <- col2rgb(w_denscols)[1,] + 1L
  w$col <- cols[w$dens]
  
  plot(phenotype ~ age, data = w[order(w$dens),], col = col,  pch = 16, 
       main = plot_title, 
       cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
       ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
       xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
  
  dd <- pdat[pdat$gender_F1M2 == 1,]
  lines(pred ~ age, data = dd, col = w_palette[8], lwd = 2)
  polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(w_palette[8], 100), border = NA)
  
  #
  # Men
  #
  m <- merged_tab[merged_tab2$gender_F1M2 == 2,]
  m_denscols <- densCols(m$age, m$phenotype, colramp=colorRampPalette(c("black", "white")))
  cols <-  colorRampAlpha(m_palette, n = 256, alpha = 0.4)
  m$dens <- col2rgb(m_denscols)[1,] + 1L
  m$col <- cols[m$dens]
  
  
  points(phenotype ~ age, data = m[order(m$dens),], col = col,  pch = 16, 
         cex = 0.6)
  
  dd <- pdat[pdat$gender_F1M2 == 2,]
  lines(pred ~ age, data = dd, col = m_palette[8], lwd = 2)
  polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(m_palette[8], 100), border = NA)
  
  
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2$age), col = "grey90")
}


draw_contour_plot <- function(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot = T, plot_title = NULL,  ymax_hist = 1, label = "", ylims_usr = NULL){
  cat("making a contour plot!\n")
  ylims <- with(merged_tab, range(phenotype))
  pheno_name <- paste(strwrap(pheno_name, width = 40), collapse = "\n")
  ylabel <- pheno_name
  pheno_name <- gsub(" \\(.*", "", pheno_name)
  binaryPhenotype <- F
  if (all(ylims == c(0,1))) {
    ylabel <- paste0(pheno_name, " frequency")
    binaryPhenotype <- T
  }
  
  merged_tab2 <- merged_tab[merged_tab$phenotype <= ylims[2] & merged_tab$phenotype >= ylims[1],]
  if (plot_title == ""){
    if (add_inter_p_to_plot) {
      if (gam.p == 0) {
        plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P < 2.23e-308")
      } else {
        plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P = ", format(gam.p, digits = 3))
      }
    } else {
      plot_title <- pheno_name
    }
  }
  if (! binaryPhenotype){
      p <- ggplot(merged_tab2) + 
        stat_density2d(geom="density2d", aes(x = age, y = phenotype, color = gender_F1M2,alpha=..level..),  contour=TRUE)+ 
        geom_line(data = pdat, aes(x = age, y = pred, color = gender_F1M2), size = 0.8) + 
        geom_ribbon(data = pdat, aes(x = age, ymin = lwr, ymax = upr, fill = gender_F1M2)) + 
        theme_minimal() +
        scale_fill_manual(values = alpha(c("#ff9999", "#99ccff"),  0.4)) + 
        theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 10), axis.text = element_text(size = 10), axis.title = element_text(size = 10)) + 
        xlab("age") + 
        ylab(ylabel) + 
        ggtitle(plot_title) + 
        scale_color_manual(values = c("indianred1", "dodgerblue1"))
  } else {
    h.br.step <- 3
    h.br <- seq(min_age, max_age, h.br.step)
    
    independ.w<- merged_tab2[merged_tab2$gender_F1M2 == 1, "age"]
    depend.w <- merged_tab2[merged_tab2$gender_F1M2 == 1, 1]
    independ.m <- merged_tab2[merged_tab2$gender_F1M2 == 2, "age"]
    depend.m <- merged_tab2[merged_tab2$gender_F1M2 == 2, 1]
    
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
    
    n_bins <- length(h.br)-1
    res_hist <- data.frame(matrix(nrow = 2*n_bins))
    res_hist$start <- rep(h.br[-length(h.br)], 2)
    res_hist$end <- rep(h.br[-1], 2)
    res_hist$gender_F1M2 <- c(rep("1", n_bins) , rep("2", n_bins))
    res_hist$val <- c(h.w, h.m)
    res_hist$bottom <- rep(0, 2*n_bins)
    
    
    p <- ggplot() + 
      geom_rect(data = res_hist, aes(ymin = bottom, ymax = val, xmin = start, xmax = end, fill = gender_F1M2, color =  gender_F1M2))+
      geom_line(data = pdat, aes(x = age, y = pred, color = gender_F1M2), size = 0.8) + 
      geom_ribbon(data = pdat, aes(x = age, ymin = lwr, ymax = upr, fill = gender_F1M2)) + 
      theme_minimal() +
      scale_fill_manual(values = alpha(c("#ff9999", "#99ccff"),  0.4)) + 
      theme(legend.position="none", plot.title = element_text(hjust = 0.5, size = 10), axis.text = element_text(size = 10), axis.title = element_text(size = 10)) + 
      xlab("age") + 
      ylab(ylabel) + 
      ylim(0,1) +
      ggtitle(plot_title) + 
      scale_color_manual(values = c("indianred1", "dodgerblue1"))
    
  }
  
  return(p)
}


# for plotting ggplots using cowplot 
format_title <- function(title_text){
  title <- ggdraw() + 
    draw_label(
      title_text,
      fontface = 'bold',
      size = 16,
      x = 0,
      hjust = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 7)
    )
  print(title)
}
