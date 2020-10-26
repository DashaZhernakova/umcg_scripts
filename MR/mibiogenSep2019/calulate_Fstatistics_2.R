args <- commandArgs(trailingOnly = TRUE)
library(TwoSampleMR)
library(MRInstruments)



get_r_from_lor <- function(dat)
{
  
  nsnp <- length(dat$SNP)
  r <- array(NA, nsnp)
  for(i in 1:nsnp)
  {
    #popaf <- get_population_allele_frequency(af[i], ncase[i] / (ncase[i] + ncontrol[i]), exp(lor[i]), prevalence[i])
    popaf <- dat$eaf.exposure
    lor = dat$beta.exposure
    
    vg <- lor[i]^2 * popaf[i] * (1-popaf[i])
    ve <- pi^2/3
    r[i] <- vg / (vg + ve)
    r[i] <- sqrt(r[i])
  }
  return(r)
}

calulateVPE <- function(exp_dat){
  VPE=(2*exp_dat$beta.exposure^2*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure))/( 2*exp_dat$beta.exposure^2*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure) + 2* exp_dat$se.exposure^2*exp_dat$samplesize.exposure*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure))
  return (VPE)
}

calculateFstats2 <- function(dat, binary = NA){
  
  if (is.na(binary)){
    if ("units.exposure" %in% colnames(dat) & !is.na(dat$units.exposure[1])){
      if (dat$units.exposure[1] == "log odds"){
        binary = T
      }
    } else if ("units.exposure_dat" %in% colnames(dat) & !is.na(dat$units.exposure_dat[1])){
      if (dat$units.exposure_dat[1] == "log odds"){
        binary = T
      }
    } else if (exp_dat$id.exposure [1] == "1112"){
      binary = T
    } else{
      print(paste("ERROR! Binary status impossible to determine!!!"))
    }
  }
  
  if (binary){
    R2 <- sum(get_r_from_lor(dat)^2)
    print(paste("Getting R2 from LOR: ",R2 ))
  } else{
    R2 <- sum(calulateVPE(dat))
    print(paste("Getting R2 from betas: ",R2 ))
  }
  print(paste("binary=", binary))
  N = mean(as.numeric(dat[, "samplesize.exposure"]))
  k = nrow(dat)
  Fstat = (R2 * (N - 1 - k))/((1 - R2) * k)
  return(Fstat)
}

getMbAF <- function(exp_dat, af_table){
  bac <- exp_dat$exposure[1]
  eafs <- rep(NA, nrow(exp_dat))
  for (i in 1:nrow(exp_dat)){
    eaf = NA
    snp <- exp_dat$SNP[i]
    snp_line <- af_table[af_table$bac == bac & af_table$rsid == snp,]
    if (snp_line$effect_allele == exp_dat$effect_allele.exposure[i] & snp_line$other_allele == exp_dat$other_allele.exposure[i]){
      eafs[i] <- snp_line$maf
    } else if (snp_line$effect_allele == exp_dat$other_allele.exposure[i] & snp_line$other_allele == exp_dat$effect_allele.exposure[i]){
      eafs[i] <- 1 - snp_line$maf
    } else{
      print(paste("Alleles not corresponding to exp_dat:", exp_dat, snpaf_table))
      eafs[i] <-  NA
    }
  }
  return(eafs)
}

get1kgAF <- function(snps, ea, oa, exp_id){
  afs <- rep(NA, length(snps))
  
  pop <- ao[ao$id == exp_id, "population"]
  if (pop == "European"){
    kg_pop <- c("CEU", "TSI", "GBR", "IBS")
  } else if (pop == "Mixed"){
    kg_pop <- c("EUR", "EAS", "SAS")
  } else{
    kg_pop <- "ALL"
  }
  
  for (i in 1:length(snps)){
    ld_res <- NULL
    tryCatch({
      ld_res <- LDpair(var1 = snps[i], var2 = snps[i], pop = kg_pop, token = "e08c86f49ed7")
    }, error=function(e) NULL)
    if (!is.null(ld_res)){
      if (ea[i] == ld_res$var1_a1 & oa[i] == ld_res$var1_a2){
        afs[i] <- ld_res$var1_a1_freq
      } else if (ea[i] == ld_res$var1_a2 & oa[i] == ld_res$var1_a1){
        afs[i] <- ld_res$var1_a2_freq
      } else{
        afs[i] <- NA
      }
    }
  }
  return (afs)
}

readExpDatFromFile <- function(exp_id, snps, exp_whole_table, af_table = NA){
  exp_table <- exp_whole_table[exp_whole_table$SNP %in% snps,]
  exp_dat <- format_data(exp_table, type = "exposure")
  #exp_dat$eaf.exposure <- getMbAF(exp_dat, af_table)
  
  return(exp_dat)
}

readExpDatFromMRbase <- function(exp_id, snps, exp_whole_table){
  exp_dat <- exp_whole_table[exp_whole_table$SNP %in% snps,]
  
  if (all(is.na(exp_dat$eaf.exposure))){
    exp_dat$eaf.exposure <- get1kgAF(exp_dat$SNP, exp_dat$effect_allele.exposure, exp_dat$other_allele.exposure, exp_dat$id.exposure[1])
  }
  #print(exp_dat$eaf.exposure)
  return(exp_dat)
}

readExpDatFromUKBFile <- function(snps, exp_whole_table){
  exp_table <- exp_whole_table[exp_whole_table$SNP %in% snps,]
  exp_dat <- format_data(exp_table, type = "exposure")
  return(exp_dat)
}
run_for_mb <- function(in_table, binary_exp = F){
  #mb_af_fpath <- args[2]
  #mb_af_fpath <- "/groups/umcg-lld/scr01/dasha/MR/results/mibiogenSep2019/fstat/snps.mafs.txt"
  #af_table <- read.table(mb_af_fpath, sep = "\t", check.names = F, as.is = T, header = T, quote="")
  
  in_table$F_exposure_2 <- NA
  all_bac <- unique(in_table$exposure)
  
  for (bac in all_bac){
    print(paste("##### Starting with", bac))
    table_subs <- in_table[in_table$exposure == bac,]
    exp_dat_path <- paste0("/groups/umcg-lld/scr01/dasha/MR/data/mibiogenOct2019/", bac, ".summary.afs.txt.gz")
    exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
    for (i in 1:nrow(table_subs)){
      outcome <- table_subs[i, "outcome"]
      print(paste(bac, outcome))
      snps <- strsplit(table_subs[i, "SNP"], "," , fixed = TRUE)[[1]]
      fstat <- NA
      tryCatch({
        exp_dat <- readExpDatFromFile(bac, snps, exp_whole_table)
        fstat <- calculateFstats2(exp_dat, binary_exp)
      }, error=function(e) {print(paste("Failed:", bac, table_subs[i, "outcome"]))})
      if (!is.na(fstat)){
        in_table[in_table$exposure == bac & in_table$outcome == outcome, "F_exposure_2"] <- fstat
        print(paste("RES", bac, outcome, fstat))
      }
    }
  }
  return (in_table)
}

run_for_gwas <- function(in_table){
  in_table$F_exposure_2 <- NA
  all_pheno_ids <- unique(in_table$id.exposure)
  
  for (pheno_id in all_pheno_ids){
    print(paste("##### Starting with", pheno_id))
    table_subs <- in_table[in_table$id.exposure == pheno_id,]
    exp_whole_table <- extract_instruments(outcomes=pheno_id, clump = FALSE)
    
    for (i in 1:nrow(table_subs)){
      outcome <- table_subs$outcome[i]
      print(paste(pheno_id, outcome))
      snps <- strsplit(table_subs[i, "SNP"], "," , fixed = TRUE)[[1]]
      fstat <- NA
      tryCatch({
        exp_dat <- readExpDatFromMRbase(pheno_id, snps, exp_whole_table)
        fstat <- calculateFstats2(exp_dat)
      }, error=function(e) {print(paste("Failed:", pheno_id, table_subs[i, "outcome"]))})
      if (!is.na(fstat)){
        in_table[in_table$id.exposure == pheno_id & in_table$outcome == outcome, "F_exposure_2"] <- fstat
      }
      print(paste("RES", pheno_id, outcome, fstat))
    }
  }
  return (in_table)
}

run_for_ukb <- function(in_table, binary){
  
  in_table$F_exposure_2 <- NA
  all_pheno <- unique(in_table$exposure)
  
  for (pheno in all_pheno){
    print(paste("##### Starting with", pheno))
    table_subs <- in_table[in_table$exposure == pheno,]
    exp_dat_path <- paste0("/groups/umcg-lld/tmp03/umcg-dzhernakova/MR/data/UKB/original_files_Neale/UKB_formatted/", pheno, ".for_MR.txt.gz")
    exp_whole_table <- read.table(gzfile(exp_dat_path), header = T, sep = "\t", as.is = T, check.names = F)
    for (i in 1:nrow(table_subs)){
      outcome <- table_subs[i, "outcome"]
      print(paste(pheno, outcome))
      snps <- strsplit(table_subs[i, "SNP"], "," , fixed = TRUE)[[1]]
      fstat <- NA
      tryCatch({
        exp_dat <- readExpDatFromUKBFile(snps, exp_whole_table)
        fstat <- calculateFstats2(exp_dat, binary)
      }, error=function(e) {print(paste("Failed:", pheno, table_subs[i, "outcome"], e))})
      if (!is.na(fstat)){
        in_table[in_table$exposure == pheno & in_table$outcome == outcome, "F_exposure_2"] <- fstat
        print(paste("RES", pheno, outcome, fstat))
      }
    }
  }
  return (in_table)
}

fpath = args[1]
exp_type = args[2]
#fpath <- "/groups/umcg-lld/scr01/dasha/MR/results/mibiogenSep2019/fstat/mb-gwas.txt"
#fpath <- "/Users/dashazhernakova/Documents/UMCG/data/MR/results2/mibiogen/mibiogenSep2019/v5/gwas-mb.v5.txt"

in_table <- read.table(fpath, sep = "\t", check.names = F, as.is = T, header = T, quote="")

if (exp_type == "mb") new_table <- run_for_mb(in_table)
if (exp_type == "gwas") new_table <- run_for_gwas(in_table)
if (exp_type == "ukb") new_table <- run_for_ukb(in_table, F)
write.table(new_table, file = paste0(fpath, ".Fstat.txt"), sep = "\t", quote = F, row.names = F)
