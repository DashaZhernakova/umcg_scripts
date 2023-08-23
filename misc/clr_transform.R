args <- commandArgs(trailingOnly = TRUE)

d <- read.delim(args[1], check.names = F, header = T, row.names = 1, as.is = T, sep ="\t")

transform_clr <- function(taxa){
  
  #Adapted from Alexander Kurilshikov 
  #replace 0 by the half of the smallest value observed
  my_min = min(taxa[taxa > 0], na.rm = T)/2
  taxa = taxa + my_min
  #Calculate geometric mean
  gm_mean = function(x, na.rm=TRUE){
    exp(mean(log(x)))
  }
  Gmean_core = apply(taxa, 1, gm_mean)
  taxa = cbind(Gmean_core, taxa)
  d <- t(apply(taxa, 1, function(b) {
    log(b/b[1])[-1]
  }))
  
  return(as.data.frame(d))
}

write.table(transform_clr(d), file = args[2], sep = "\t", quote = F, col.names = NA)