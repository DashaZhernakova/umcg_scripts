library(corrplot)
args <- commandArgs(trailingOnly = TRUE)

wd_path <- "/groups/umcg-lifelines/tmp01/users/umcg-dzhernakova/gender_difs/"

# Phenotypes
#traits_path <- "Laboratory_assessment_Blood_1A.dat"
traits_path <- args[1]
st_col = 2

pdfpath <- paste0("plots/", traits_path, "_corrplot2.pdf")
out_path <- paste0("summary_tables/", traits_path, "_correl2.txt")

setwd(wd_path)

# traits of interest
traits0 <- read.delim(traits_path, header = T, row.names = 1, sep = "\t", as.is = T, check.names = F)
traits0 <- traits0[,seq(st_col,ncol(traits0))]
traits <- sapply(traits0, function(x) as.numeric(as.character(x)))
row.names(traits) <- row.names(traits0)


cor.mtest <- function(mat, ...) {
    mat <- as.data.frame(mat)
    n <- ncol(mat)
    p.mat<- matrix(NA, n, n)
    r.mat<- matrix(NA, n, n)
    diag(p.mat) <- 0
    for (i in 1:(n - 1)) {
        for (j in (i + 1):n) {
		if ( nrow(mat[!is.na(mat[,i]) & !is.na(mat[,j]),]) < 10) {
			p.mat[i, j] <- p.mat[j, i] <- 1
			r.mat[i, j] <- r.mat[j, i] <- 0
		} else {	
		#print (paste0(colnames(mat)[i], " ", colnames(mat)[j], " ", length(mat[!is.na(mat[,i]),i]) , " ", length(mat[!is.na(mat[,j]),j]) ))
            	tmp <- cor.test(mat[, i], mat[, j], ...)
            	p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
		r.mat[i, j] <- r.mat[j, i] <- tmp$estimate
        	}
	}
    }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  colnames(r.mat) <- rownames(r.mat) <- colnames(mat)
  return(list(r.mat, p.mat))
}


res_lst <- cor.mtest(traits, method = "spearman", use = "pairwise.complete.obs")

c <- res_lst[[1]]
diag(c) <- 1
pvs <- res_lst[[2]]
#c <- cor(traits, method = "spearman", use = "pairwise.complete.obs")
write.table(c, file = out_path, sep = "\t", quote = F, col.names = NA)

pdf(pdfpath, width = 15, height = 15)
corrplot(c, type="upper",  p.mat = pvs, sig.level = 0.05, order = 'hclust', insig = "blank")
dev.off()
