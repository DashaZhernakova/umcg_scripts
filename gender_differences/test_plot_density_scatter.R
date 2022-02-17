colorRampAlpha <- function(..., n, alpha) {
  colors <- colorRampPalette(...)(n)
  paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
}




draw_smooth_scatter <- function(merged_tab, pheno_name, pdat, gam.p, min_age, max_age, add_inter_p_to_plot = T, plot_title = NULL, plot_points = T, breakpoints = NULL, factor_name = "", alpha_points = 40, breakpoints_intervals = NULL, ymax_hist = 1, label = "", ylims_usr = NULL){
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
  
  plot(phenotype ~ age, data = merged_tab2,  col = gender_F1M2,  pch = 16, 
       main = plot_title, 
       cex = 0.6, xlab = "age", ylab = ylabel, cex.main = cex_main, frame.plot = F, axes = T, 
       ylim =c(min(pretty(merged_tab2$phenotype)), max(pretty(merged_tab2$phenotype))),
       xlim = c(min(min_age,merged_tab2$age), max(max_age, merged_tab2$age)))
  
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
  polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(w_palette[8], 50), border = NA)
  
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
  
  dd <- pdat[pdat$gender_F1M2 == 1,]
  lines(pred ~ age, data = dd, col = m_palette[8], lwd = 2)
  polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(w_palette[8], 50), border = NA)
  
  
  abline(h = pretty(merged_tab2$phenotype), col = "grey90")
  abline(v = pretty(merged_tab2$age), col = "grey90")
}