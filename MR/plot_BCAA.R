args <- commandArgs(trailingOnly = TRUE)
library(TwoSampleMR)
library(cowplot)
library(ggplot2)

#
# !!! Modify the paths !!!
#

access_token_fname = "C:/Users/Dasha/work/UMCG/data/MR/results2/mrbase.oauth"
plot_dir = "C:/Users/Dasha/work/UMCG/data/MR/results2/AA_T2D/plots/"
setwd("C:/Users/Dasha/work/UMCG/data/MR/results2/AA_T2D/")
pval_thres = 5e-08
#
#

custom_scatter_plot <- function(res, dat){
  d <- subset(dat, mr_keep)
  if(nrow(d) < 2 | sum(d$mr_keep) == 0)
  {
    return(blank_plot("Insufficient number of SNPs"))
  }
  
  index <- d$beta.exposure < 0
  d$beta.exposure[index] <- d$beta.exposure[index] * -1
  d$beta.outcome[index] <- d$beta.outcome[index] * -1
  mrres <- subset(res, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
  mrres$a <- 0
  exp_name <- gsub(' \\|\\|.*',"",d$exposure[1])
  out_name <- gsub(' \\|\\|.*',"",d$outcome[1])
  
  text_x <- max(d$beta.exposure + d$se.exposure)
  text_y <- min(d$beta.outcome-d$se.outcome)
  
  cols = c("red1", "dodgerblue1")
  p <- ggplot(data=d, aes(x=beta.exposure, y=beta.outcome)) +
    geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
    geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
    geom_point(aes(text=paste("SNP:", SNP)), colour = "grey30") +
    geom_abline(data=mrres, aes(intercept=a, slope=b), colour = cols[2], size = 1, show.legend=F) +
    labs(colour="MR Test", x=paste("SNP effect on", exp_name), y=paste("SNP effect on", out_name)) +
    theme_bw() +
    annotate("text", text_x, text_y, size = 3, label = paste0(format(mrres$b, digits = 2), " SD change in ", exp_name, "\nper SD change in ", out_name), hjust = 1, vjust = 0)
  return(p)
}

run_mr_and_plot <- function(exp, out){
  exp_from_file = is.na(as.numeric(exp))
  out_from_file = is.na(as.numeric(out))
  
  print(paste("Exposure = ", exp))
  print(paste("Outcome = ", out))
  
  print(paste("Exposure will be read from file = ", exp_from_file))
  print(paste("Outcome will be read from file = ", out_from_file))
  
  
  if (exp_from_file & !out_from_file){ 
    # if exposure should be read from file and outcome - from MRBase
    exp_dat_path <- exp
    outcome_id <- out
    
    exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
    exp_table <- exp_whole_table[exp_whole_table$pval < as.numeric(pval_thres),]
    out_dat <- extract_outcome_data(snps = exp_table$SNP, outcomes = outcome_id, access_token = access_token_fname)
    exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
  } else if (!exp_from_file & out_from_file){
    # if exposure is in MRBase, and outcome - in a file
    exp_id <- exp
    out_dat_path <- out
    out_whole_table <- read.table(gzfile(out_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
    
    exp_dat <- extract_instruments(outcomes=exp_id, clump = FALSE, access_token = access_token_fname)
    out_dat <- format_data(out_whole_table, snps = exp_dat$SNP, type = "outcome")
    exp_dat <- clump_data(exp_dat[exp_dat$SNP %in% out_dat$SNP,])
  } else if (exp_from_file & out_from_file){
    #if both are in files
    exp_dat_path <- exp
    out_dat_path <- out
    
    exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
    exp_table <- exp_whole_table[exp_whole_table$pval < as.numeric(pval_thres),]
    
    out_whole_table <- read.table(gzfile(out_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
    
    out_dat <- format_data(out_whole_table, snps = exp_table$SNP, type = "outcome")
    exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
  } else{
    # if both are in MRBase
    exp_id <- exp
    outcome_id <- out
    
    exp_dat <- extract_instruments(outcomes = exp_id, access_token = access_token_fname)
    out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id, access_token = access_token_fname)
  }
  
  dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
  res <- mr(dat, method_list=c("mr_ivw", "mr_wald_ratio"))
  
  p <- custom_scatter_plot(res, dat)
  return(p)
}


exp_from_file = FALSE
out_from_file = FALSE

plot_fname = paste0(plot_dir, "selected_MR_plots_v3.pdf")

to_plot <- matrix(c("2", "897", "49", "897", "2", "940", "936", "873", "850", "873", "940", "igfbp2.txt.gz"), ncol = 2, byrow = T)
plots <- list()

for (i in 1:nrow(to_plot)){
  p <- run_mr_and_plot(to_plot[i,1], to_plot[i,2])
  plots[[i]] <- p
}


p <- plot_grid(plotlist = plots, nrow = 2,  axis = 'lb', align = 'hv', scale = 0.9)
save_plot(plot_fname, p, ncol = 3, nrow = 2, base_asp = 1.2)

print("Finished plotting!")