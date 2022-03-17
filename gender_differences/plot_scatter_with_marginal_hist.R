add_prevalence_hist_to_plot <- function(merged_tab, scale.hist = 5, las.h1 = 1, ylabel2 = "Count"){
  
  col.hist.w = col2transparent("indianred1",50)
  col.hist.m = col2transparent("dodgerblue1", 50)
  
  label = paste0("Cohen's f^2 = ", formatC(gam.cohen_f2, digits = 4, format = "f"))
  pheno_name <- paste(strwrap(pheno_name, width = 40), collapse = "\n")
  ylabel <- pheno_name
  pheno_name <- gsub(" \\(.*", "", pheno_name)

 
  max_pheno <- max(merged_tab$phenotype)
  min_pheno <- min(merged_tab$phenotype)
  
  n_breaks <- 20
  h.br.x.step <- 3
  h.br.y.step <- (max_pheno - min_pheno)/n_breaks
  h.br.x <- seq(min_age, max_age, by = h.br.x.step)
  h.br.y <- seq(min_pheno, max_pheno, by = h.br.y.step)
  
  age.w<- merged_tab[merged_tab$gender_F1M2 == 1, "age"]
  pheno.w <- merged_tab[merged_tab$gender_F1M2 == 1, 1]
  age.m <- merged_tab[merged_tab$gender_F1M2 == 2, "age"]
  pheno.m <- merged_tab[merged_tab$gender_F1M2 == 2, 1]
  
  # X axis (age)
  h.x <- hist(merged_tab$age,
              breaks = h.br.x,
              plot = FALSE
  )$mid

  #men - age
  h.m.x <- hist(age.m,
                 breaks = h.br.x,
                 plot = FALSE
  )$counts

  # women - age
  h.w.x <- hist(age.w,
                 breaks = h.br.x,
                 plot = FALSE
  )$counts
  
  # Y axis (pheno)
  h.y <- hist(merged_tab$phenotype,
              breaks = h.br.y,
              plot = FALSE
  )$mid
  #men
  h.m.y <- hist(pheno.m,
                breaks = h.br.y,
                plot = FALSE
  )$counts
  
  # women
  h.w.y <- hist(pheno.w,
                breaks = h.br.y,
                plot = FALSE
  )$counts
  
  
  ymin <- min(merged_tab$phenotype)
  ymax <- max(merged_tab$phenotype)
  xmin <- min(min_age, merged_tab$age)
  xmax <- max(max_age, merged_tab$age)
  if (gam.p == 0) {
    plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P < 2.23e-308")
  } else {
    plot_title <- paste0(pheno_name, '\n', label, "\nGAM interaction P = ", format(gam.p, digits = 3))
  }
  pdf(paste0(plot_path, ".test.pdf"))
  #def.par <- par(no.readonly = TRUE)
  par(mar = c(1,1,1,1))
  #mar = c(1.5,2,4,1),
  layout(matrix(c(1,1,1,3,2,2,3,2,2,5,4,4), 4, 3, byrow = TRUE))
  
  plot.new()
  plot.window(0:1, 0:1)
  text(0.5, 0.5, plot_title, cex = 3)
  
  par(mgp = c(1.2, 0.4, 0),
      tck = -.01, # Reduce tick length
      xaxs = "i", yaxs = "i", las = 1) # Remove plot padding
  palette(c(col2transparent("#ff9999", 120),col2transparent("#99ccff", 120)))
  plot(phenotype ~ age, data = merged_tab,  col = gender_F1M2,  pch = 16, 
       xlab = NA, ylab = NA,
       cex = 0.6,  frame.plot = F, axes = T, 
       ylim =c(ymin, ymax),
       xlim = c(xmin, xmax))
  abline(h = pretty(merged_tab$phenotype), col = "grey90")
  abline(h = min(merged_tab$phenotype), col = "grey90")
  abline(v = pretty(merged_tab$age), col = "grey90")
  
  levs <- levels(merged_tab$gender_F1M2)
  cols = c("indianred1", "dodgerblue1")
  
  ## add the fitted lines
  for (l in seq_along(levs)) {
    dd <- pdat[pdat$gender_F1M2 == levs[l],]
    lines(pred ~ age, data = dd, col = cols[[l]], lwd = 2)
    polygon(c(rev(dd$age), dd$age), c(rev(dd$lwr), dd$upr), col = col2transparent(cols[[l]], 100), border = NA)
  }
  
  
  # draw histogram:
  
  #barplot(-h.m.y, xaxt='n', space = 0, horiz = T, ylim = c(ymin, ymax), border = "darkgrey", col = col2transparent("dodgerblue1", 25))
  #barplot(-h.w.y, xaxt='n', horiz = T, border = "darkgrey", col = col2transparent("indianred1", 25), add = T)

  
  # phenotype hist
  par(mar = c(1,1,1,1.6))
  plot (NA, type='n', axes = F, xaxt='n', 
        xlab=NA, ylab=NA, main=NA,
        ylim=c(min_age,max_age),
        xlim=c(-max(h.m.y, h.w.y), 0))
  mtext(ylabel, 4, las = 3, cex = 0.8)
  
  for (i in 1:length(h.y)) {
    if (h.m.y[i] > 0) {
      polygon(c(0, rep(-1*h.m.y[i], 2), 0),
              c(rep(h.br.x[i], 2), rep(h.br.x[i + 1], 2)),
              col = col.hist.m, border = "darkgrey"
      )
    }
    if (h.w.y[i] > 0) {
      polygon(c(0, rep(-1*h.w.y[i], 2), 0),
              c(rep(h.br.x[i], 2), rep(h.br.x[i + 1], 2)),
              col = col.hist.w, border = "darkgrey"
      )
    }
  }
  
  # age hist
  par(mar = c(1,1,1.6,1))
  plot (NA, type='n', axes = F, xaxt='n',
        xlab=NA, ylab=NA, main=NA,
        xlim=c(min_age,max_age),
        ylim=c(-max(h.m.x, h.w.x), 0))
  mtext("age", 3, cex = 0.8)
  
  for (i in 1:length(h.x)) {
    if (h.m.x[i] > 0) {
      polygon(c(rep(h.br.x[i], 2), rep(h.br.x[i + 1], 2)),
              c(0, rep(-h.m.x[i], 2), 0),
              col = col.hist.m, border = "darkgrey"
      )
    }
    if (h.w.x[i] > 0) {
      polygon(c(rep(h.br.x[i], 2), rep(h.br.x[i + 1], 2)),
              c(0, rep(-h.w.x[i], 2), 0),
              col = col.hist.w, border = "darkgrey"
      )
    }
  }
  
  
}
  