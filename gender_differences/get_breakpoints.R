library(ggpmisc)

get_peak_locations_pval <- function(x, span_val = 7, pval_threshold = 5e-4){
  pow <- floor(log10(pval_threshold))
  pval_thres2 <- 5*10^(pow/2) 
  peaks <- ggpmisc:::find_peaks(-x, span = span_val) & x < pval_threshold & c(x[2:length(x)], 1) < pval_thres2 & c(1, x[1:(length(x) - 1)]) < pval_thres2
  return(peaks)
}

get_peak_locations <- function(x, span_val = 7, tval_threshold = 4){
  tval_thres2 <- 2
  x <- abs(x)
  peaks <- ggpmisc:::find_peaks(x, span = span_val) & x > tval_threshold & c(x[2:length(x)], 1) > tval_thres2 & c(1, x[1:(length(x) - 1)]) > tval_thres2
  return(peaks)
}

get_breakpoints <- function(merged_tab, ttest_window = 5){
  #merged_tab <- merged_tab[merged_tab$age <= 80 & merged_tab$age > 20,]
  #age2test <- (min(merged_tab$age) + window):(max(merged_tab$age) - window)
  age2test = (max(20, min(merged_tab$age)) + ttest_window) : (min(68, max(merged_tab$age)) - ttest_window)
  res <- data.frame(matrix(ncol = 2, nrow = length(age2test)))
  colnames(res) <- c("w", "m")
  rownames(res) <- as.character(age2test)
  #p_thres <- ifelse (nrow(merged_tab) > 5000, 1e-30, 5e-4) 
  p_thres <- 4
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
  
  peaks_m_idx <- rownames(res)[get_peak_locations(res[,"m"], span_val = 7)]
  peaks_w_idx <- rownames(res)[get_peak_locations(res[,"w"], span_val = 7)]
  
  return (list(peaks_w_idx, peaks_m_idx))
}
