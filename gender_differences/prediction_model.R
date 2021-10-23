library(mgcv)

setwd("/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/prediction_model")
d <- read.delim("age_gender_bmi_glu_smk_lipids_sbp_t2d_antihyper_newcvd_all_LL_with_children.txt", header = T, row.names = 1, as.is = T, check.names = F, sep = "\t")
d <- na.omit(d)
m1 <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) + SMK3 + s(CHO) + s(HDC) + s(SBP), data = d, method = "ML", family = binomial(link = "logit"))

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

m2 <- gam(new_cvd~gender_F1M2 + s(age) + s(BMI) + SMK3 + s(CHO) + s(HDC) + s(SBP) + ti(BMI, age, by = gender_F1M2) + ti(CHO, age, by =gender_F1M2) + ti(HDC, age, by = gender_F1M2) + ti(SBP, age, by =gender_F1M2), data = d, method = "ML", family = binomial(link = "logit"))
summary(m2)
#Parametric coefficients:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept) -4.53682    0.06627 -68.455  < 2e-16 ***
#gender_F1M2  0.27054    0.03999   6.765 1.34e-11 ***
#SMK3         0.19512    0.04494   4.341 1.42e-05 ***
#Approximate significance of smooth terms:
#                          edf Ref.df   Chi.sq  p-value
#s(age)                  3.756  4.660 1431.115  < 2e-16 ***
#s(BMI)                  1.000  1.000    7.997 0.004692 **
#s(CHO)                  4.259  5.188   53.107 6.04e-10 ***
#s(HDC)                  2.124  2.706   23.188 2.68e-05 ***
#s(SBP)                  1.000  1.000   36.558 1.49e-09 ***
#ti(BMI,age):gender_F1M2 1.001  1.002    0.000 0.984363
#ti(CHO,age):gender_F1M2 1.982  2.438   30.258 9.47e-07 ***
#ti(HDC,age):gender_F1M2 2.400  2.912    9.131 0.024103 *
#ti(SBP,age):gender_F1M2 2.340  2.841   20.333 0.000107 ***


summary(m1)
#Parametric coefficients:
#            Estimate Std. Error z value Pr(>|z|)
#(Intercept) -4.60862    0.06490 -71.011  < 2e-16 ***
#gender_F1M2  0.28657    0.03919   7.312 2.64e-13 ***
#SMK3         0.20357    0.04484   4.540 5.63e-06 ***
#
#Approximate significance of smooth terms:
#         edf Ref.df  Chi.sq  p-value
#s(age) 4.318  5.287 1723.65  < 2e-16 ***
#s(BMI) 1.000  1.000   17.61 2.71e-05 ***
#s(CHO) 4.227  5.146   99.60  < 2e-16 ***
#s(HDC) 2.494  3.160   75.30 7.90e-16 ***
#s(SBP) 1.000  1.000   61.63 4.16e-15 ***

AIC(m1,m2)
#         df      AIC
#m1 17.27538 29286.36
#m2 26.74855 29237.87

anova