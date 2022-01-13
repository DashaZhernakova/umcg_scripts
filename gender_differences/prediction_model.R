library(mgcv)

setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/prediction_model")

traits_path <- "../v4/data/LL_phenotypes_merged_all.log_some.v5.txt"
pheno_path <- "age_gender_new_CVD.v5.txt"

traits <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
pheno <- read.table(pheno_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits_m <- traits[match(row.names(pheno), row.names(traits), nomatch = 0 ), , drop = F]
pheno_m <- pheno[match(row.names(traits_m), row.names(pheno), nomatch = 0), ]
all(row.names(traits_m) == row.names(pheno_m))
num_traits <- ncol(traits_m)
d <- cbind(pheno_m, traits_m)

m0 <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D + s(CHO) + s(HDC) + s(SBP), data = d, method = "ML", family = binomial(link = "logit"), select = T)
summary(m0)
m1 <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D + s(CHO) + s(HDC) + s(SBP) + s(LLDS_T1A) + s(total_scor_VAL), data = d, method = "ML", family = binomial(link = "logit"), select = T)
summary(m1)
m2 <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D+ s(CHO) + s(HDC) + s(SBP)  + s(LLDS_T1A) + s(total_scor_VAL) + s(age, by = gender_F1M2) + 
            ti(BMI, age, by = gender_F1M2) + ti(CHO, age, by =gender_F1M2) + ti(HDC, age, by = gender_F1M2) + ti(SBP, age, by =gender_F1M2) + ti(LLDS_T1A, age, by = gender_F1M2) + ti(total_scor_VAL, age, by =gender_F1M2), data = d, method = "ML", family = binomial(link = "logit"), select=T)

m3 <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D + s(CHO) + s(HDC) + s(SBP)  + s(LLDS_T1A) + s(total_scor_VAL) + 
            s(age, by = gender_F1M2) + s(BMI, by = gender_F1M2) + interaction(SMK3, gender_F1M2) + interaction(T2D, gender_F1M2) + s(CHO, by = gender_F1M2) + s(HDC, by = gender_F1M2) + s(SBP, by = gender_F1M2)  + s(LLDS_T1A, by = gender_F1M2) + s(total_scor_VAL, by = gender_F1M2) +
            ti(BMI, age) + s(age, by = SMK3) + s(age, by = T2D) + ti(CHO, age) + ti(HDC, age) + ti(SBP, age) + ti(LLDS_T1A, age) + ti(total_scor_VAL, age) + 
            s(age, by = interaction(SMK3, gender_F1M2)) + s(age, by = interaction(T2D, gender_F1M2)) + ti(BMI, age, by = gender_F1M2) + ti(CHO, age, by =gender_F1M2) + ti(HDC, age, by = gender_F1M2) + ti(SBP, age, by =gender_F1M2) + ti(LLDS_T1A, age, by = gender_F1M2) + ti(total_scor_VAL, age, by =gender_F1M2), data = d, method = "ML", family = binomial(link = "logit"), select=T)
print(summary(m1))
print(summary(m2))
print(summary(m3))
print(AIC(m0, m1, m2, m3))

# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept)                       -4.39189    0.11143 -39.414  < 2e-16 ***
#   gender_F1M2                        0.00000    0.00000      NA       NA
# SMK3                               0.39946    0.09742   4.101 4.12e-05 ***
#   T2D                                0.35341    0.18477   1.913   0.0558 .
# interaction(SMK3, gender_F1M2)1.1  0.00000    0.00000      NA       NA
# interaction(SMK3, gender_F1M2)0.2 -0.17499    0.26707  -0.655   0.5123
# interaction(SMK3, gender_F1M2)1.2 -0.22241    0.28821  -0.772   0.4403
# interaction(T2D, gender_F1M2)1.1   0.00000    0.00000      NA       NA
# interaction(T2D, gender_F1M2)0.2   0.13060    0.25869   0.505   0.6137
# interaction(T2D, gender_F1M2)1.2   0.00000    0.00000      NA       NA
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value
# s(age)                                   1.2046850      9 277.203  < 2e-16 ***
#   s(BMI)                                   1.1539722      9   4.490 0.019149 *
#   s(CHO)                                   1.8671186      9  11.733 0.000689 ***
#   s(HDC)                                   0.8080999      9   3.749 0.017674 *
#   s(SBP)                                   2.6085895      9  58.899  < 2e-16 ***
#   s(LLDS_T1A)                              0.0013966      9   0.001 0.360370
# s(total_scor_VAL)                        0.6757301      9   1.993 0.019129 *
#   s(age):gender_F1M2                       0.0002252      9   0.000 0.252428
# s(BMI):gender_F1M2                       0.0003890      9   0.000 1.000000
# s(CHO):gender_F1M2                       0.0007567      9   0.000 0.388593
# s(HDC):gender_F1M2                       0.0003850      9   0.000 0.645612
# s(SBP):gender_F1M2                       0.0004743      9   0.000 0.649864
# s(LLDS_T1A):gender_F1M2                  0.5297943      9   1.188 0.132145
# s(total_scor_VAL):gender_F1M2            1.3305963      9   5.083 0.005212 **
#   ti(BMI,age)                              0.0010709     16   0.000 0.654202
# s(age):SMK3                              0.0003212      9   0.000 0.318236
# s(age):T2D                               0.0001275      9   0.000 0.307392
# ti(CHO,age)                              0.0009806     16   0.001 0.532281
# ti(HDC,age)                              2.8949273     16   9.100 0.003832 **
#   ti(SBP,age)                              1.2581371     16   3.778 0.041129 *
#   ti(LLDS_T1A,age)                         0.0006983     16   0.000 0.590131
# ti(total_scor_VAL,age)                   0.0009046     16   0.000 0.745605
# s(age):interaction(SMK3, gender_F1M2)0.1 0.0005523      9   0.000 0.664255
# s(age):interaction(SMK3, gender_F1M2)1.1 0.0005419      9   0.000 0.668238
# s(age):interaction(SMK3, gender_F1M2)0.2 0.0002697      9   0.000 0.514987
# s(age):interaction(SMK3, gender_F1M2)1.2 0.0006295      9   0.000 0.384618
# s(age):interaction(T2D, gender_F1M2)0.1  0.0005010      9   0.000 0.858928
# s(age):interaction(T2D, gender_F1M2)1.1  0.0001893      9   0.000 0.762190
# s(age):interaction(T2D, gender_F1M2)0.2  1.7647189      9   6.233 0.014372 *
#   s(age):interaction(T2D, gender_F1M2)1.2  0.0001622      9   0.000 0.316803
# ti(BMI,age):gender_F1M2                  0.0009178     16   0.000 0.595443
# ti(CHO,age):gender_F1M2                  0.0019982     16   0.002 0.412901
# ti(HDC,age):gender_F1M2                  1.5405915     16   3.691 0.021934 *
#   ti(SBP,age):gender_F1M2                  0.0009576     16   0.000 0.915181
# ti(LLDS_T1A,age):gender_F1M2             0.0007453     16   0.000 0.789820
# ti(total_scor_VAL,age):gender_F1M2       1.1121596     16   2.016 0.122607
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# Rank: 423/427
# R-sq.(adj) =  0.025   Deviance explained = 9.09%
# -ML = 6983.5  Scale est. = 1         n = 69159
# df      AIC
# m0 18.44275 16643.48
# m1 22.30745 13979.04
# m2 25.45681 13971.71
# m3 30.94226 13961.28




summary(m0)
# Family: binomial
# Link function: logit
# 
# Formula:
#   new_cvd ~ gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D + s(CHO) +
#   s(HDC) + s(SBP)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -4.47797    0.08628 -51.903  < 2e-16 ***
#   gender_F1M2  0.17340    0.05271   3.290   0.0010 **
#   SMK3         0.32282    0.05942   5.433 5.56e-08 ***
#   T2D          0.21182    0.12233   1.732   0.0834 .
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value
# s(age) 4.7266      9 713.000  < 2e-16 ***
#   s(BMI) 0.8588      9   6.080  0.00752 **
#   s(CHO) 1.5812      9   6.291  0.01608 *
#   s(HDC) 1.7769      9  26.747 7.62e-08 ***
#   s(SBP) 2.6164      9  86.243  < 2e-16 ***
#   ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# R-sq.(adj) =  0.0232   Deviance explained = 8.41%
# -ML = 8325.5  Scale est. = 1         n = 82345

summary(m1)
# Family: binomial
# Link function: logit
# 
# Formula:
#   new_cvd ~ gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D + s(CHO) +
#   s(HDC) + s(SBP) + s(LLDS_T1A) + s(total_scor_VAL)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -4.49818    0.09459 -47.552  < 2e-16 ***
#   gender_F1M2  0.17078    0.05791   2.949  0.00319 **
#   SMK3         0.36532    0.06533   5.592 2.24e-08 ***
#   T2D          0.20942    0.13239   1.582  0.11369
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq p-value
# s(age)            4.4851      9 596.597 < 2e-16 ***
#   s(BMI)            1.1613      9   4.936 0.01655 *
#   s(CHO)            1.8535      9  10.352 0.00188 **
#   s(HDC)            2.1227      9  25.946 1.9e-07 ***
#   s(SBP)            2.5167      9  79.367 < 2e-16 ***
#   s(LLDS_T1A)       0.3043      9   0.427 0.23030
# s(total_scor_VAL) 1.7038      9  15.333 6.1e-05 ***
#   ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# R-sq.(adj) =  0.0244   Deviance explained = 8.87%
# -ML = 6992.7  Scale est. = 1         n = 69159


summary(m2)
# 
# Family: binomial
# Link function: logit
# 
# Formula:
#   new_cvd ~ gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D + s(CHO) +
#   s(HDC) + s(SBP) + s(LLDS_T1A) + s(total_scor_VAL) + s(age,
#                                                         by = gender_F1M2) + ti(BMI, age, by = gender_F1M2) + ti(CHO,
#                                                                                                                 age, by = gender_F1M2) + ti(HDC, age, by = gender_F1M2) +
#   ti(SBP, age, by = gender_F1M2) + ti(LLDS_T1A, age, by = gender_F1M2) +
#   ti(total_scor_VAL, age, by = gender_F1M2)
# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -4.47644    0.09517 -47.037  < 2e-16 ***
#   gender_F1M2  0.15636    0.05878   2.660  0.00782 **
#   SMK3         0.36959    0.06535   5.656 1.55e-08 ***
#   T2D          0.20128    0.13288   1.515  0.12983
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value
# s(age)                             2.1629273      9 519.432  < 2e-16 ***
#   s(BMI)                             1.1756710      9   5.132  0.01446 *
#   s(CHO)                             1.7393476      9   7.706  0.00432 **
#   s(HDC)                             1.7008735      9   8.612  0.00238 **
#   s(SBP)                             2.4770921      9  56.359  < 2e-16 ***
#   s(LLDS_T1A)                        0.3479688      9   0.529  0.21614
# s(total_scor_VAL)                  1.7452812      9  13.716 6.74e-05 ***
#   s(age):gender_F1M2                 0.0001033      9   0.000  0.08135 .
# ti(BMI,age):gender_F1M2            0.0007760     16   0.000  0.83168
# ti(CHO,age):gender_F1M2            0.6650100     16   0.888  0.23415
# ti(HDC,age):gender_F1M2            1.7047875     16   7.926  0.00493 **
#   ti(SBP,age):gender_F1M2            0.6562781     16   0.937  0.20924
# ti(LLDS_T1A,age):gender_F1M2       0.0013165     16   0.000  0.72934
# ti(total_scor_VAL,age):gender_F1M2 1.3269443     16   2.239  0.14331
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# R-sq.(adj) =  0.0245   Deviance explained = 8.95%
# -ML = 6986.6  Scale est. = 1         n = 69159





# old summary(m2)

# 
# Parametric coefficients:
#   Estimate Std. Error z value Pr(>|z|)
# (Intercept) -4.45957    0.08632 -51.665  < 2e-16 ***
#   gender_F1M2  0.16606    0.05296   3.135  0.00172 **
#   SMK3         0.32737    0.05936   5.515 3.49e-08 ***
#   T2D          0.21197    0.12265   1.728  0.08395 .
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# Approximate significance of smooth terms:
#                           edf Ref.df  Chi.sq  p-value
#   s(age)                  2.4857095      9 646.269  < 2e-16 ***
#   s(BMI)                  0.8735221      9   6.881  0.00457 **
#   s(CHO)                  1.3748965      9   3.780  0.04404 *
#   s(HDC)                  1.2709529      9   4.297  0.02673 *
#   s(SBP)                  2.6256786      9  84.787  < 2e-16 ***
#   s(age):gender_F1M2      0.0007349      9   0.003  0.04025 *
#   ti(BMI,age):gender_F1M2 0.0020910     16   0.000  1.00000
#   ti(CHO,age):gender_F1M2 0.7211479     16   1.083  0.18825
#   ti(HDC,age):gender_F1M2 2.0627389     16  18.119 1.91e-05 ***
#   ti(SBP,age):gender_F1M2 0.0083449     16   0.010  0.26656
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# R-sq.(adj) =  0.0231   Deviance explained = 8.49%
# -ML = 8316.2  Scale est. = 1         n = 82345



summary(m1)
#Parametric coefficients:
# Estimate Std. Error z value Pr(>|z|)
# (Intercept) -4.47797    0.08628 -51.903  < 2e-16 ***
#   gender_F1M2  0.17340    0.05271   3.290   0.0010 **
#   SMK3         0.32282    0.05942   5.433 5.56e-08 ***
#   T2D          0.21182    0.12233   1.732   0.0834 .
# ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df  Chi.sq  p-value
#   s(age) 4.7266      9 713.000  < 2e-16 ***
#   s(BMI) 0.8588      9   6.080  0.00752 **
#   s(CHO) 1.5812      9   6.291  0.01608 *
#   s(HDC) 1.7769      9  26.747 7.62e-08 ***
#   s(SBP) 2.6164      9  86.243  < 2e-16 ***
#   ---
#   Signif. codes:  0 ---***--- 0.001 ---**--- 0.01 ---*--- 0.05 ---.--- 0.1 --- --- 1
# 
# R-sq.(adj) =  0.0232   Deviance explained = 8.41%
# -ML = 8325.5  Scale est. = 1         n = 82345

m0 <- gam(new_cvd~gender_F1M2 + age + BMI + SMK3 + T2D + CHO + HDC + SBP, data = d, method = "ML", family = binomial(link = "logit"), select = T)


AIC(m1,m2)
#df      AIC
#m1 18.44275 16643.48
#m2 18.82591 16628.35


anova(m2,m1, test="Chisq")
# P = 5.415e-05

preds <- colnames(d)[c(4:14, 16:48,50,51)]
terms <-  paste0(" + s(", paste(preds, collapse = ")+ s("), ")")
m3 <- gam(as.forumla(paste("new_cvd~gender_F1M2 + s(age) + s(BMI) + SMK3 + T2D+ s(age, by=gender_F1M2)", terms, sep = " ")), data = d, method = "ML", family = binomial(link = "logit"), select=T) 







m_bmi <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) , data = d, method = "ML", family = binomial(link = "logit"))
m_bmi_full <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) +ti(BMI, age, by = gender_F1M2) , data = d, method = "ML", family = binomial(link = "logit"))
AIC(m_bmi, m_bmi_full)

m_cho <- gam(new_cvd~gender_F1M2 + s(age) + s(CHO) , data = d, method = "ML", family = binomial(link = "logit"))
m_cho_full <- gam(new_cvd~gender_F1M2 + s(age) + s(CHO) +ti(CHO, age, by =gender_F1M2) , data = d, method = "ML", family = binomial(link = "logit"))
AIC(m_cho, m_cho_full)

m_hdc <- gam(new_cvd~gender_F1M2 + s(age) + s(HDC) , data = d, method = "ML", family = binomial(link = "logit"))
m_hdc_full <- gam(new_cvd~gender_F1M2 + s(age) + s(HDC) +ti(HDC, age, by = gender_F1M2) , data = d, method = "ML", family = binomial(link = "logit"))
summary(m_hdc_full)
AIC(m_hdc, m_hdc_full)

m_sbp <- gam(new_cvd~gender_F1M2 + s(age) + s(SBP) , data = d, method = "ML", family = binomial(link = "logit"))
m_sbp_full <- gam(new_cvd~gender_F1M2 + s(age) + s(SBP) +ti(SBP, age, by =gender_F1M2) , data = d, method = "ML", family = binomial(link = "logit"))
summary(m_sbp_full)
AIC(m_sbp, m_sbp_full)

m_smk <- gam(new_cvd~gender_F1M2 + s(age) + SMK3 , data = d, method = "ML", family = binomial(link = "logit"))
m_smk_full <- gam(new_cvd~gender_F1M2 + s(age) + SMK3 + s(age, by = interaction(SMK3, gender_F1M2)) , data = d, method = "ML", family = binomial(link = "logit"))
summary(m_smk_full)
AIC(m_smk, m_smk_full)

