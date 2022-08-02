args <- commandArgs(trailingOnly = TRUE)
library(TwoSampleMR)
library(MRInstruments)
library("LDlinkR")

ao <- available_outcomes()

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
  #VPE=(2*exp_dat$beta.exposure^2*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure))/( 2*exp_dat$beta.exposure^2*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure) + 2* exp_dat$se.exposure^2*exp_dat$samplesize.exposure*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure))
  VPE=(2*exp_dat$beta.exposure^2*exp_dat$eaf.exposure*(1-exp_dat$eaf.exposure))
  
  return (VPE)
}

calculateFstats2 <- function(dat, binary = NA){
  
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

get1kgAF <- function(snps, ea, oa, exp_id){
  afs <- rep(NA, length(snps))
  
  #pop <- ao[ao$id == exp_id, "population"]
  pop = "European"
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