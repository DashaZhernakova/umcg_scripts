x <- read.delim("../all_LL/no_menopause/f2_cmp.txt", sep = "\t",  as.is = T, header = 1)

pdf("../all_LL/no_menopause/f2_cmp.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
max_sex_f2 <- max(x$sex_gam_f2_before, x$sex_gam_f2_after)
plot(x$sex_gam_f2_before, x$sex_gam_f2_after, pch = 16, 
     xlim = c(0, max_sex_f2), ylim = c(0, max_sex_f2),
     ylab = "Cohen's f2 of sex for individuals above 55 years old", xlab = "Cohen's f2 of sex for individuals below 45 years old",
     main = "Comparison of the effect size of sex\nbefore and after menopause")
abline(c(0,0), c(1,1), lty = 2)
subs <-x[x$sex_gam_f2_before > 0.01 | x$sex_gam_f2_after > 0.01,]
points(subs$sex_gam_f2_before, subs$sex_gam_f2_after, pch = 16, 
       xlim = c(0, max_sex_f2), ylim = c, col = "red")

max_inter = max(x$cohen_f2_before, x$cohen_f2_after)
plot(x$cohen_f2_before, x$cohen_f2_after, pch = 16, 
     xlim = c(0, max_inter), ylim = c(0, max_inter),
     ylab = "Cohen's f2 of age by sex interaction after 55", xlab = "Cohen's f2 of age by sex interaction before 45",
     main = "Comparison of the age by sex interaction\neffect size before and after menopause")
abline(c(0,0), c(1,1), lty = 2)

max_age = max(x$age_gam_f2_before, x$age_gam_f2_after)
plot(x$age_gam_f2_before, x$age_gam_f2_after, pch = 16, 
     xlim = c(0, max_age), ylim = c(0, max_age),
     ylab = "Cohen's f2 of age after 55", xlab = "Cohen's f2 of age before 45",
     main = "Comparison of the effect size of age\nbefore and after menopause")
abline(c(0,0), c(1,1), lty = 2)
dev.off()

