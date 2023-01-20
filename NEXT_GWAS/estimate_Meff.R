setwd("/Users/Dasha/work/UMCG/data/NEXT/metabolites/")


fname <- "metabo_Baby_B.filtered.txt"
d <- as.data.frame(t(read.delim(fname, header = T, sep = "\t", row.names = 1, as.is = T, check.names = F)))
corr.matrix<-abs(cor(d))

evals<-eigen(t(corr.matrix),symmetric=T)$values
oldV<-var(evals)
M<-length(evals)
L<-(M-1)
Meffold<-M*(1-(L*oldV/M^2))

if (evals == 1) { 
  oldV <- 0 
  Meffold <- M
}

labelevals<-array(dim=M)
for(col in 1:M) { labelevals[col]<-c(col) }
levals<-cbind(labelevals, evals)

newevals<-evals
for(i in 1:length(newevals)) { 
  if(newevals[i] < 0) { 
    newevals[i] <- 0
  }
}

newlevals<-cbind(labelevals, newevals)

newV<-var(newevals)
Meffnew<-M*(1-(L*newV/M^2))
Meffnew
