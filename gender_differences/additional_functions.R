library(RColorBrewer)
library('dplyr')

calculate_sex_diff_ttest <- function(merged_tab, covariates = c(), min_age = 20, max_age = 80){
  colnames(merged_tab)[1] <- "phenotype"
  merged_tab <- merged_tab[(merged_tab$age < max_age) & (merged_tab$age >= min_age),]
  merged_tab <- mutate(merged_tab, gender_F1M2 = factor(gender_F1M2))
  
  if (length(covariates) == 0){
    resid <- merged_tab$phenotype 
    #res <- lm(phenotype ~ gender_F1M2, data = merged_tab)
  } else{
    resid <- residuals( lm(as.formula(paste("phenotype ~ ",  paste(covariates, collapse = "+"))), data = merged_tab))
  }
  #sex_p <- summary(res)$coefficients[2,4]
  t_res <- t.test(resid[merged_tab$gender_F1M2 == 1], resid[merged_tab$gender_F1M2 == 2])
  sex_p <- t_res$p.value
  return(sex_p)
}

test_polynomial_interaction <- function(merged_tab, model_fam = NO(), gm_mean_cov = F){
  if (! gm_mean_cov){
    #lgam.fit <- gamlss(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) + I(age^4), data = merged_tab, family = model_fam, method = RS(30))
    #lm2 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^3) + I(age^4) + gender_F1M2:(age + I(age^2) + I(age^3) + I(age^4)), data = merged_tab, family = model_fam, method = RS(30)) 
    lgam.fit <- gamlss(phenotype ~ gender_F1M2 + age + I(age^2) + I(age^3) , data = merged_tab, family = model_fam, method = RS(30))
    lm2 <- gamlss(phenotype ~ gender_F1M2 + age + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)), data = merged_tab, family = model_fam, method = RS(30)) 
  } else {
    lgam.fit <- gamlss(phenotype ~ gender_F1M2 + age + gm_mean + I(age^2) + I(age^3) , data = merged_tab, family = model_fam, method = RS(30))
    lm2 <- gamlss(phenotype ~ gender_F1M2 + age + gm_mean + I(age^3) + gender_F1M2:(age + I(age^2) + I(age^3)), data = merged_tab, family = model_fam, method = RS(30)) 
    
  }
  
  lr <- LR.test(lgam.fit, lm2, print = F)
  m <- lgam.fit
  if (lr$p.val < 0.05){
    m <- lm2
  }
  return (list("lr_p.val" = lr$p.val, "mod" = m))
}


# Calculates Cohen's R^2 to estimate the local effect size of the variable of interest
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3328081/
# According to Cohen�s (1988) guidelines, f2>=0.02, f2>=0.15, and f2>=0.35 represent small, medium, and large effect sizes, respectively.
calculate_cohens_f2 <- function(mod_full, mod_part){
  r_full <- summary(mod_full)$r.sq
  r_sub <- summary(mod_part)$r.sq
  f2 <- (r_full - r_sub)/(1-r_full)
  return(f2)
}

correct_for_covariates_before <- function(traits_m, pheno_m, covariates){
  new_traits <- matrix(nrow = nrow(traits_m), ncol = ncol(traits_m))
  colnames(new_traits) <- colnames(traits_m)
  rownames(new_traits) <- rownames(traits_m)
  
  if (! all(covariates %in% colnames(pheno_m))){
    print("Error! Some covariates are not present in the phenotype file")
    print(paste0("Covariates: ", paste(covariates, collapse = ",")))
    print(paste0("Covariate file phenotypes: ", paste(colnames(pheno_m), collapse = ",")))
  }
  
  for (col in 1:ncol(traits_m)){
    d <- cbind(traits_m[,col], pheno_m[,covariates])
    colnames(d)[1] <- "phenotype"
    lm.fit <- lm(as.formula(paste("phenotype ~  ", paste(covariates, collapse = "+"), sep = " ")), data = d)
    new_traits[,col] <- residuals(lm.fit)
  }
  return(new_traits)
}


run_for_split_by_covariate <- function(merged_tab, pheno_name, covariate_to_split, highlight_positive = F, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots = T,  gam_family = gaussian(), min_age = 20, max_age = 80){
  pdat_list <- list()
  col_list <- list()
  colors <- list("indianred1", "orange1", "dodgerblue1", "cadetblue1")
  w0 <- merged_tab[merged_tab$gender_F1M2 == 1 & merged_tab[,covariate_to_split] == 0,]
  w1 <- merged_tab[merged_tab$gender_F1M2 == 1 & merged_tab[,covariate_to_split] == 1,]
  m0 <- merged_tab[merged_tab$gender_F1M2 == 2 & merged_tab[,covariate_to_split] == 0,]
  m1 <- merged_tab[merged_tab$gender_F1M2 == 2 & merged_tab[,covariate_to_split] == 1,]  
  
  cnt <- 1
  cnt2 <- 1
  for (dat in list(w0,w1,m0,m1)){
    if (nrow(dat) > 0){
      pdat_list[[cnt]] <- fit_gam_age_only(dat, pheno_name, covariates_linear = c(), covariates_nonlinear = c(), n_points = 300, make_plots = T,  gam_family = gaussian(), min_age = 20, max_age = 80)
      col_list[[cnt]] <- colors[[cnt2]]
      cnt <- cnt + 1
    }
    cnt2 <- cnt2 + 1
  }
  factor_to_highlight <- ""
  if (highlight_positive) factor_to_highlight <- covariate_to_split
  colnames(merged_tab)[1] <- "phenotype"
  draw_plot_multiline(merged_tab, pheno_name, pdat_list, col_list, factor_name = factor_to_highlight)
}
