#
# Various phenotypes
#


setwd("/Users/dashazhernakova/Documents/UMCG/data/gender_differences/omics/anova+gam")

# NMR metabolomics
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/Metabolomics/metabolomics_shared_folder/2.metabolites/LLD_bloodlipids_nmr.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- traits[,seq(5,ncol(traits))]

# proteins
traits <- as.data.frame(t(read.table("/Users/dashazhernakova/Documents/UMCG/data/olink/olink_shared_folder/Data/rawProteinData/CVD3_olinkNormal_1447_LLDsamples_t_ProtNames.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))
# drugs
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/LifeLines_phenotypes/20150618_45drugs_1135metasubjects_atleast5users.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
# pheno from Science
traits <- read.table("C:/Users/Dasha/work/UMCG/data/LifeLines_phenotypes/20150715_intristic_1135patients.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- traits[,seq(3,ncol(traits))]
# pheno more samples
traits <- read.table("/Users/dashazhernakova/Documents/UMCG/data/LifeLines_phenotypes/20170123_selection_phenotypes_for_TL_not_binary.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits <- traits[,seq(3,ncol(traits))]

cvd <- c("relativesCVD","hr","p","p_axis","pq","qrs","qrs_axis","qt","qtc","t_axis","dbp","hbf","map","sbp","cho","crph","glu","hb1c","hbac","hdc","k","ldc","tgl","ins","homair","bmi","hip","angioplasty","arythmia","chestpain","diabetes","diabtype","dilataorta","edema","heartattack","heartfailure","highchol","hypertension","lof","narrowcarotis","stroke","T1Dy_n","T2Dy_n","RelativesDiab","added_compl_legs_1_yes_0_no","added_pain_hands_or_feet_1y_0n","added_stiff_hands_or_feet_1y_0n","feetwounds","jointpainfeet","jointpainhands","legcomplnight","legcomplwalking","stiffnessfeet","stiffnesshands")
#traits <- as.data.frame(t(read.table(gzfile("/Users/dashazhernakova/Documents/UMCG/data/olink/LLDeep_expression/gene_level/LLD_genelevel_readcounts.nozeros.TMM.Log2Transformed.ProbesCentered.SamplesZTransformed.txt.gz"), header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

num_traits <- ncol(traits)

pheno0 <- as.data.frame(t(read.table("age_gender_smk_contrac_cell_counts.cc_outliers_na.txt", header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)))

pheno <- na.omit(pheno0)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ),]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0),]
all(row.names(traits_m) == row.names(pheno_m))

traits_m <- traits_m[order(pheno_m$age),]
pheno_m <- pheno_m[order(pheno_m$age),]



##########################


make_plots = T
pdf("../anova+gam/rnaseq_plots/ageDEgenes/new/phenotypes.pdf", width = 15, height = 15)
par(mfrow=c(5,4)) 

n_points = 300
res_dif_all2 <- data.frame(age = seq(20, 75, length = n_points))
res_summary <- data.frame()

cnt = 1

indices = 1:ncol(traits_m)
#indices = sample(1:ncol(traits_m), 500, replace=F)
#plot_list <- list()

for (idx in indices){
  
  print(idx)
  
  trait_id <- colnames(traits_m)[idx]
  trait_name = trait_id
  if (length(trait_name) > 0){ #if gene id in gene convertion table
    merged_tab <- rm_na_outliers(traits_m, pheno_m, idx)
    
    res_dif = NULL
    res_dif_lst <- plot_scatter_and_gam2(merged_tab, trait_name, T, 300, make_plots)
    
    if (res_dif_lst[["inter_p"]] < 0.05){
      cnt <- cnt + 1
      res_dif_all2[,trait_id] <- res_dif_lst[["dif"]]
      
      res_summary[trait_id,'inter_p'] = res_dif_lst[["inter_p"]]
      res_summary[trait_id,'g_beta'] = res_dif_lst[["g_beta"]]
      res_summary[trait_id,'g_pv'] = res_dif_lst[["g_pv"]]
      
      #res_summary[trait_id,'p1'] = res_dif_lst[["plots"]][[1]]
      #res_summary[trait_id,'p2'] = res_dif_lst[["p2"]]
      #plot_list[[idx]] = res_dif_lst[["plots"]]
    }
  }
  
}

pheno_diff_all <- as.data.frame(t(res_dif_all2))
pheno_name = "Biochem_Cholesterol"
for (g in row.names(t_res_dif_all)){
  c <- cor(as.numeric(t_res_dif_all[g,]), as.numeric(pheno_diff_all[pheno_name,]))
  g_name = gene_table[gene_table[,1] == g,2]
  if(c  > 0.95){
    print(paste(g, g_name, c))
  }
}




