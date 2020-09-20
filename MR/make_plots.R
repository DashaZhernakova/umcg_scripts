args <- commandArgs(trailingOnly = TRUE)
library(TwoSampleMR)
library(cowplot)
library(ggplot2)

#
# !!! Modify the paths !!!
#

access_token_fname = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/mrbase.oauth"
plot_dir = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/mibiogenSep2019/plots/"

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
  
  cols = c("red1", "dodgerblue1")
  ggplot(data=d, aes(x=beta.exposure, y=beta.outcome)) +
    geom_errorbar(aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
    geom_errorbarh(aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
    geom_point(aes(text=paste("SNP:", SNP))) +
    geom_abline(data=mrres, aes(intercept=a, slope=b), colour = ifelse(mrres$b > 0, cols[1], cols[2]), size = 1.2, show.legend=F) +
    labs(colour="MR Test", x=paste("SNP effect on", exp_name), y=paste("SNP effect on", out_name)) +
    theme_bw() +
    annotate("text", Inf, Inf, label = paste0(format(mrres$b, digits = 3), " SD change in ", exp_name, "\nper SD change in ", out_name), hjust = 1, vjust = 1)
  
}

exp_from_file = FALSE
out_from_file = FALSE

plot_fname = paste0(plot_dir,gsub(".*/", "", args[1]), "_vs_", gsub(".*/", "", args[2]), ".pdf")

exp_from_file = is.na(as.numeric(args[1]))
out_from_file = is.na(as.numeric(args[2]))

print(paste("Exposure = ", args[1]))
print(paste("Outcome = ", args[2]))
print(paste("P-value threshold for exposure = ", args[3]))
print(paste("Output plot path = ", plot_fname))

print(paste("Exposure will be read from file = ", exp_from_file))
print(paste("Outcome will be read from file = ", out_from_file))


if (exp_from_file & !out_from_file){ 
  # if exposure should be read from file and outcome - from MRBase
  exp_dat_path <- args[1]
  outcome_id <- args[2]
  pval_thres = args[3]
  
  exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  exp_table <- exp_whole_table[exp_whole_table$pval < as.numeric(pval_thres),]
  out_dat <- extract_outcome_data(snps = exp_table$SNP, outcomes = outcome_id, access_token = access_token_fname)
  exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
  } else if (!exp_from_file & out_from_file){
  # if exposure is in MRBase, and outcome - in a file
  exp_id <- args[1]
  out_dat_path <- args[2]
  out_whole_table <- read.table(gzfile(out_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  
  exp_dat <- extract_instruments(outcomes=exp_id, clump = FALSE, access_token = access_token_fname)
  out_dat <- format_data(out_whole_table, snps = exp_dat$SNP, type = "outcome")
  exp_dat <- clump_data(exp_dat[exp_dat$SNP %in% out_dat$SNP,])
} else if (exp_from_file & out_from_file){
  #if both are in files
  exp_dat_path <- args[1]
  out_dat_path <- args[2]
  pval_thres = args[3]
  
  exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  exp_table <- exp_whole_table[exp_whole_table$pval < as.numeric(pval_thres),]
  
  out_whole_table <- read.table(gzfile(out_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  
  out_dat <- format_data(out_whole_table, snps = exp_table$SNP, type = "outcome")
  exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
} else{
  # if both are in MRBase
  exp_id <- args[1]
  outcome_id <- args[2]
  
  exp_dat <- extract_instruments(outcomes = exp_id, access_token = access_token_fname)
  out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id, access_token = access_token_fname)
}

dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)
res <- mr(dat, method_list=c("mr_ivw", "mr_wald_ratio"))

plots = list()

p1 <- mr_scatter_plot(res, dat)
plots[[1]] = p1[[1]]

res_single <- mr_singlesnp(dat, all_method = c("mr_ivw", "mr_weighted_median", "mr_egger_regression"))
p2 <- mr_forest_plot(res_single)
plots[[2]] = p2[[1]]

res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
plots[[3]] = p3[[1]]

p4 <- mr_funnel_plot(res_single)
plots[[4]] = p4[[1]]

p <- plot_grid(plotlist = plots, nrow = 2,  axis = 'lb', align = 'hv', labels=c('Scatter plot', 'Forest plot', 'Leave-one-out plot', 'Funnel plot'), scale = 0.9)
save_plot(plot_fname, p, ncol = 2, nrow = 2)

print("Finished plotting!")