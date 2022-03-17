covar <- read.delim("data/LLD.covariates.txt",header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
sex <- read.delim("../genotypes/LLD/LLD_filtered.fam", header = F, sep = "\t", row.names = 2)
covar <- covar[,c("age", "F.prausnitzii", "read_number")]

svs <- read.delim("data/LLD.dSV.filtered.txt",header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
clades <- read.delim("data/Fprau_clades.txt",header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
ids <- intersect(row.names(svs), row.names(clades))
sv_names <- c("F.prausnitzii:102", "F.prausnitzii:101", "F.prausnitzii:69", "F.prausnitzii:9")

svs[] <- lapply(as.data.frame(svs), function(x) sub(1,0,x))
svs[] <- lapply(as.data.frame(svs), function(x) sub(2,1,x))
for (sv in sv_names){
    for (cl in colnames(clades)){
        #m <- as.data.frame(cbind(svs[ids,sv], clades[ids,cl]))
        #m <- na.omit(m)
        #cat(sv, cl, chisq.test(table(m))$p.value, "\n", sep = " ")
        m <- as.data.frame(cbind(svs[ids,sv], clades[ids,cl], covar[ids,], sex[ids, "V5"]))
        colnames(m) <- c("sv", "clade", "age", "abundance", "read_number", "sex")
        pval <- summary(glm(sv ~ clade + age + sex + abundance + read_number, data = m, family = binomial(link = "logit")))$coefficients["clade",4]
        cat(sv, cl, pval, "\n", sep = "\t")
    }
}

svs <- read.delim("data/LLD.vSV.filtered.txt",header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
clades <- read.delim("data/Fprau_clades.txt",header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
ids <- intersect(row.names(svs), row.names(clades))

sv="F.prausnitzii:33"
for (cl in colnames(clades)){
    #m <- as.data.frame(cbind(svs[ids,sv], clades[ids,cl]))
    #m <- na.omit(m)
    #cat(sv, cl, wilcox.test(m[,1] ~ m[,2])$p.value, "\n", sep = " ")
    m <- as.data.frame(cbind(svs[ids,sv], clades[ids,cl], covar[ids,], sex[ids, "V5"]))
    colnames(m) <- c("sv", "clade", "age", "abundance", "read_number", "sex")
    pval <- summary(lm(sv ~ clade + age + sex + abundance + read_number, data = m))$coefficients["clade",4]
    cat(sv, cl, pval, "\n", sep = "\t")
}

