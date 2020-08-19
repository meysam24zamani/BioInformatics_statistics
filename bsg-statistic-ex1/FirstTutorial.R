load(url("http://www-eio.upc.es/~jan/data/bsg/Chromosome1_CHBPopSubset.rda"))

dim(Ysub)
Ysub[1:5, 1:5]

install.packages("genetics")
library(genetics)

Ysub[Ysub == "NN"] = NA

snp1 <- Ysub[,1]
snp1 <- genotype(snp1, sep="")
summary(snp1)

snp2 <- Ysub[,2]
snp2 <- genotype(snp2, sep="")
summary(snp2)

snp3 <- Ysub[,3]
snp3 <- genotype(snp3, sep="")
summary(snp3)

pmis <- function(x) {
  n <- length(x)
  nmis <- sum(is.na(x))
  y <- 100*nmis/n
  return(y)
}

pmis(Ysub[,1])
pmis(Ysub[,2])

per.mis.persnp <- apply(Ysub, 2, pmis)

plot(1:ncol(Ysub), per.mis.persnp, ylim = c(0,100))

per.mis.perind <- apply(Ysub, 1, pmis)

plot(1:nrow(Ysub), per.mis.perind, ylim = c(0,100))

n <- nrow(Ysub)
p <- ncol(Ysub)
sum(is.na(Ysub))/(n*p)*100

summary(snp3)

snp3
table(snp3)
nCC <- sum(snp3=="C/C", na.rm = T)
nTC <- sum(snp3=="C/T" | snp3=="T/C", na.rm = T)
nTT <- sum(snp3=="T/T", na.rm = T)
nCC
nTC
nTT

nC <- 2*nCC + nTC
nT <- 2*nTT + nTC
nC
nT

n
nmis <- sum(is.na(snp3))
nTotal <- 2*n - 2*nmis

nTotal

nC+nT

pC <- nC/nTotal
pT <- nT/nTotal

pC+pT
pC
pT
res <- summary(snp3)

min(pC, pT)

attributes(res)

res$allele.freq[,2]

maf <- function(x) {
  res <- summary(x)
  af <- res$allele.freq[,2]
  minaf <- min(af, na.rm=T)
  
  if (minaf==1) 
    minaf <- 0
  
  return(minaf)
}

y <- apply(Ysub, 2, maf)
