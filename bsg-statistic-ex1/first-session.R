rm(list = ls())
ls()
load(url("Chromosome1_CHBPopSubset.rda")

dim(Ysub) Ysub[1:5,5:1]    

install.packages("genetics")
library(genetics)

Ysub[Ysub="NN"] <- NA

snpl <- Ysub[,1]
snp1 <- genotype(snp1,sep="")
summary(snp1)

snp2 <- Ysub[,2]
snp2 <- genotype(snp2,sep="")
summary(snp2)

snp3 <- Ysub[,3]
snp3 <- genotype(snp3,sep="")
summary(snp3)

pmis <- function(x){
  n <- length(x)
  nmis <- sum(is.na(x))
  y <- 100*nmis/n
  return(y)
}

pmis(Ysub[,1])
pmis(Ysub[,2])

per.mis.persnp <- apply(Ysub,2,pmis)
plot(1:ncol(Ysub),per.mis.persnp,ylim=c())




snap3
table(snap3)
nCC <- sum(snp3=="C/C")
nTC <- sum(snp3=="C/T" | snp3=="T/C",a.rm=TRUE)
nTT <- sum(snp3=="T/T",na.rm=TRUE)
nCC


PC <- NC/nTotal
res <- summary(np3)
min(pC,pT)


maf <- function(x){

  
}