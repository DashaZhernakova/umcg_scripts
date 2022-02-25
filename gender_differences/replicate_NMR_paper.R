indices = 1:ncol(traits_m)
cnt = 1
cat("\nStarting the analyses\n")
for (idx in indices){
  
  trait_name <- ifelse(is.null(pheno_table), colnames(traits_m)[idx], pheno_table[pheno_table[,1] == colnames(traits_m)[idx], 2])
  
  log_transform = FALSE
  
  merged_tab <- rm_na_outliers(traits_m, pheno_m, idx, method = outlier_correction_method, log_tr = log_transform, scale_tr = scale_transform)
  res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, covariates_linear = covariateslinear2, covariates_nonlinear = covariatesnonlinear2, n_points = n_points, make_plots = F, gam_family = gam_family, min_age = min_age, max_age = max_age, ymax_hist = ymax_hist, label = '', add_inter_p_to_plot = add_inter_p_to_plot, plot_title = plot_title, interp_cutoff = interp_cutoff, plot_points = plot_points, add_breakpoints = add_breakpoints,  t_threshold = ttest_cutoff, derivatives_cutoff = deriv_cutoff, log_tr = log_transform)
  pdat <- res_dif_lst$pdat
  
  res25 <- 0
  res50 <- 0
  if (! is.null(pdat)){
    pdat50m <- pdat[pdat$age > 49 & pdat$age < 51 & pdat$gender_F1M2 == 2,]
    pdat50w <- pdat[pdat$age > 49 & pdat$age < 51 & pdat$gender_F1M2 == 1,]
    
    if (mean(pdat50m$pred) > mean(pdat50w$pred)){
      if (min(pdat50m$lwr) > max(pdat50w$upr)){
        res50 <- 1
      }
    } else if (mean(pdat50w$pred) > mean(pdat50m$pred)){
      if (min(pdat50w$lwr) > max(pdat50m$upr)){
        res50 <- -1
      }
    }
    pdat25m <- pdat[pdat$age > 24 & pdat$age < 26 & pdat$gender_F1M2 == 2,]
    pdat25w <- pdat[pdat$age > 24 & pdat$age < 26 & pdat$gender_F1M2 == 1,]
    
    if (mean(pdat25m$pred) > mean(pdat25w$pred)){
      if (min(pdat25m$lwr) > max(pdat25w$upr)){
        res25 <- 1
      }
    } else if (mean(pdat50w$pred) > mean(pdat50m$pred)){
      if (min(pdat50w$lwr) > max(pdat50m$upr)){
        res25 <- -1
      }
    }
  }
  cat(idx, colnames(traits_m)[idx], trait_name, res25, res50, "\n", sep = "\t") 
}
  