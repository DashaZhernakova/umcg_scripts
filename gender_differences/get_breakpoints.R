library(ggpmisc)

get_peak_locations_pval <- function(x, span_val = 7, pval_threshold = 5e-4){
  pow <- floor(log10(pval_threshold))
  pval_thres2 <- 5*10^(pow/2) 
  peaks <- ggpmisc:::find_peaks(-x, span = span_val) & x < pval_threshold & c(x[2:length(x)], 1) < pval_thres2 & c(1, x[1:(length(x) - 1)]) < pval_thres2
  return(peaks)
}

get_peak_locations <- function(x, span_val = 4, tval_threshold = 4){
  tval_thres2 <- 2
  x <- abs(x)
  peaks <- ggpmisc:::find_peaks(x, span = span_val) & x > tval_threshold & c(x[2:length(x)], 1) > tval_thres2 & c(1, x[1:(length(x) - 1)]) > tval_thres2
  return(peaks)
}

get_breakpoints_ttest <- function(merged_tab, ttest_window = 5, tval_threshold = 4){
  #merged_tab <- merged_tab[merged_tab$age <= 80 & merged_tab$age > 20,]
  age2test <- (max(min(merged_tab[merged_tab$gender_F1M2 == 1,"age"]), min(merged_tab[merged_tab$gender_F1M2 == 2,"age"])) + 3):
    (min(max(merged_tab[merged_tab$gender_F1M2 == 1,"age"]), max(merged_tab[merged_tab$gender_F1M2 == 2,"age"])) - 3)
  #age2test = (max(20, min(merged_tab$age)) + ttest_window) : (min(68, max(merged_tab$age)) - ttest_window)
  res <- data.frame(matrix(ncol = 2, nrow = length(age2test)))
  colnames(res) <- c("w", "m")
  rownames(res) <- as.character(age2test)

  for (i in age2test){
    left <- merged_tab[merged_tab$age <= i & merged_tab$age > (i - ttest_window),]
    right <- merged_tab[merged_tab$age > i & merged_tab$age <= (i + ttest_window),]

    t_w <- t.test(left[left$gender_F1M2 == 1,1], right[right$gender_F1M2 == 1,1])
    t_m <- t.test(left[left$gender_F1M2 == 2,1], right[right$gender_F1M2 == 2,1])
    res[as.character(i),"w"] <- t_w$statistic
    res[as.character(i),"m"] <- t_m$statistic
    #print(paste(i, t_w$statistic, t_w$p.value, t_m$statistic, t_m$p.value))
  }
  #res
  
  peaks_m_idx <- rownames(res)[get_peak_locations(res[,"m"], span_val = 7, tval_threshold = tval_threshold)]
  peaks_w_idx <- rownames(res)[get_peak_locations(res[,"w"], span_val = 7, tval_threshold = tval_threshold)]
  
  return (list(peaks_w_idx, peaks_m_idx))
}

get_breakpoints_within_interval <- function(deriv_breakpoints, ttest_breakpoints, age_interval1 = 10, age_interval2 = 3){
  res <- list()
  visited = c()
  cnt <- 1
  for (der in as.numeric(deriv_breakpoints)) {
    for (br in as.numeric(ttest_breakpoints)) {
      #if (! br %in% visited){
        if ((br - der < age_interval1) & (br > der)){
          res[[cnt]] <- c(der, br)
          cnt <- cnt + 1
          visited <- c(visited,br)
        } else if ((der - br < age_interval2) & (der > br)){
          res[[cnt]] <- c(br, der)
          cnt <- cnt + 1
          visited <- c(visited,br)
        }
      #}
    }
  }
  return (res)
}

get_breakpoints <- function(merged_tab, correctForCellCounts, t_threshold = 4, derivatives_cutoff = 0.0003, min_age = 20, max_age = 80, age_interval1 = 10, age_interval2 = 10){
  ttest_breakpoints <- get_breakpoints_ttest(merged_tab, tval_threshold = t_threshold,ttest_window = 10)
  deriv_breakpoints <- get_breakpoints_derivatives(merged_tab, correctForCellCounts = correctForCellCounts, cutoff = derivatives_cutoff, min_age = min_age, max_age = max_age)
  
  res_w <- get_breakpoints_within_interval(deriv_breakpoints[[1]], ttest_breakpoints[[1]], age_interval1 = age_interval1, age_interval2 = age_interval2)
  res_m <- get_breakpoints_within_interval(deriv_breakpoints[[2]], ttest_breakpoints[[2]], age_interval1 = age_interval1, age_interval2 = age_interval2)
  
  return(list(res_w, res_m))
}
