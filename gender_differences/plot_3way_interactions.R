library(gratia)

colnames(merged_tab)[1] <- "phenotype"
w <- merged_tab[merged_tab$gender_F1M2 == 1,]
m <- merged_tab[merged_tab$gender_F1M2 == 2,]

gam_fit <- gam(phenotype~s(age) + s(ly) + s(bmi) + ti(age,bmi) , 
               data = merged_tab, method = "REML")
summary(gam_fit)

gam_fit_m <- gam(phenotype~s(age) + s(ly) + s(bmi) + ti(age,bmi) , 
               data = m, method = "REML")
summary(gam_fit_m)

gam_fit_w <- gam(phenotype~s(age) + s(ly) + s(bmi) + ti(age,bmi) , 
               data = w, method = "REML")
summary(gam_fit_w)


w <- merged_tab[merged_tab$gender_F1M2 == 1,]
m <- merged_tab[merged_tab$gender_F1M2 == 2,]
gam_fit_w <- gam(phenotype ~   s(age) + s(stress_1y) + s(stress_chronic) + s(phys_activity_total) + s(phys_activity_intensive) + s(diet) + s(alcohol) + ti(diet, age), method = "REML", select = T, data = w)
summary(gam_fit_w)
draw(gam_fit_w)
gam_fit_m <- gam(phenotype ~   s(age) + s(stress_1y) + s(stress_chronic) + s(phys_activity_total) + s(phys_activity_intensive) + s(diet) + s(alcohol) + ti(diet, age), method = "REML", select = T, data = m)
gam_fit <- gam(phenotype ~   s(age) + s(stress_1y) + s(stress_chronic) + s(phys_activity_total) +  s(phys_activity_intensive) + s(diet) + s(alcohol) + ti(diet, age), data = merged_tab, method = "REML", select = T)
summary(gam_fit)
summary(ga)

