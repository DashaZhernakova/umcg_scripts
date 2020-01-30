args <- commandArgs(trailingOnly = TRUE)

library(TwoSampleMR)
library(MRInstruments)

#
# !!! Modify the paths !!!
#
source("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/run_MR.R")
access_token_fname = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/v2/mrbase.oauth"
out_dir = "/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/results/mibiogenSep2019/per_bact3/"

#
#

exp_from_file = FALSE
out_from_file = FALSE

exp_from_file = is.na(as.numeric(args[1]))
out_from_file = is.na(as.numeric(args[2]))

print(paste("Exposure = ", args[1]))
print(paste("Outcome = ", args[2]))
print(paste("P-value threshold for exposure = ", args[3]))

print(paste("Exposure will be read from file = ", exp_from_file))
print(paste("Outcome will be read from file = ", out_from_file))


if (exp_from_file & !out_from_file){ 
  # if exposure should be read from file and outcome - from MRBase
  exp_dat_path <- args[1]
  outcome_id <- args[2]
  pval_thres = args[3]
  
  exp_name <- basename(exp_dat_path)
  out_name <- outcome_id
  
  exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  exp_table <- exp_whole_table[exp_whole_table$pval < as.numeric(pval_thres),]
  out_dat <- extract_outcome_data(snps = exp_table$SNP, outcomes = outcome_id, access_token = access_token_fname)
  exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
  
  } else if (!exp_from_file & out_from_file){
  # if exposure is in MRBase, and outcome - in a file
  exp_id <- args[1]
  out_dat_path <- args[2]
  pval_thres = 5e-08
  out_whole_table <- read.table(gzfile(out_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  
  exp_name <- exp_id
  out_name <- basename(out_dat_path)
  
  exp_dat <- extract_instruments(outcomes=exp_id, clump = FALSE, access_token = access_token_fname)
  out_dat <- format_data(out_whole_table, snps = exp_dat$SNP, type = "outcome")
  exp_dat <- clump_data(exp_dat[exp_dat$SNP %in% out_dat$SNP,])
  
} else if (exp_from_file & out_from_file){
  #if both are in files
  exp_dat_path <- args[1]
  out_dat_path <- args[2]
  pval_thres = args[3]
  
  exp_name <- basename(exp_dat_path)
  out_name <- basename(out_dat_path)
  
  exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  exp_table <- exp_whole_table[exp_whole_table$pval < as.numeric(pval_thres),]
  
  out_whole_table <- read.table(gzfile(out_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
  
  out_dat <- format_data(out_whole_table, snps = exp_table$SNP, type = "outcome")
  exp_dat <- clump_data(format_data(exp_table, snps = out_dat$SNP, type = "exposure"))
} else{
  # if both are in MRBase
  exp_id <- args[1]
  outcome_id <- args[2]
  pval_thres = 5e-08
  
  exp_name <- exp_id
  out_name <- out_id
  
  exp_dat <- extract_instruments(outcomes = exp_id, access_token = access_token_fname)
  out_dat <- extract_outcome_data(snps = exp_dat$SNP, outcomes = outcome_id, access_token = access_token_fname)
}

out_path <- paste0(out_dir, pval_thres, "/")
out_fname <- paste0(out_path, exp_name, "-", out_name, ".txt")

res <- run_mr(exp_dat, out_dat)

if (!is.null(res)){
  write.table(res, file =  out_fname , sep = "\t", quote = F, col.names = NA)
}
print("Finished")