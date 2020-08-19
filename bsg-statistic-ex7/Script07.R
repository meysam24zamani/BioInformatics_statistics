#
# Script session 7: relatedness investigation
#

rm(list=ls())

library(data.table)


X <- fread("http://www-eio.upc.es/~jan/Data/bsg/YRI06.raw",data.table = FALSE)

class(X)
X[1:5,1:10]
X.Fam <- X[,1:6]
X.Gen <- X[,7:ncol(X)]
X <- as.matrix(X.Gen)

ibs.mean <- function(x,y) {
   y <- mean(2 - abs(x - y),na.rm=TRUE)
   return(y)
}

ibs.sd <- function(x,y) {
   y <- sd(abs(x-y),na.rm=TRUE)
   return(y)
}

ibs.mean(X[1,],X[2,])
ibs.sd(X[1,],X[2,])

ibs.mean(X[1,],X[1,])
ibs.sd(X[1,],X[1,])

n <- nrow(X)
n

p <- ncol(X)
p

Dmean <- matrix(NA,nrow=n,ncol=n)
Dsd <- matrix(NA,nrow=n,ncol=n)

#
# calculate m,s for each pair of individuals
#

for(i in 1:n) {
   for(j in 1:n) {
      Dmean[i,j] <- ibs.mean(X[i,],X[j,])
      Dsd[i,j] <- ibs.sd(X[i,],X[j,])
   }
}

Dmean[1:5,1:5]
Dsd[1:5,1:5]

ibs.m <- Dmean[lower.tri(Dmean)]
ibs.s <- Dsd[lower.tri(Dsd)]

plot(ibs.m,ibs.s,xlim=c(1,2),ylim=c(0,1),xlab="Mean",ylab="Standard deviation")

head(X.Fam)

Relationship <- matrix(NA,nrow=n,ncol=n)

is.po <- function(x,y) {
   xchildy  <- x[1]==y[1] & (x[3]==y[2] | x[4]==y[2]) # x child of y
   xparenty <- x[1]==y[1] & (x[2]==y[3] | x[2]==y[4]) # x parent of y
   y <- xchildy | xparenty
   return(y)
}

head(X.Fam)
is.po(X.Fam[1,],X.Fam[2,])
is.po(X.Fam[1,],X.Fam[3,])
is.po(X.Fam[2,],X.Fam[3,])
is.po(X.Fam[2,],X.Fam[1,])
is.po(X.Fam[3,],X.Fam[1,])
is.po(X.Fam[3,],X.Fam[2,])

for(i in 1:n) {
   for(j in 1:n) {
      Relationship[i,j] <- is.po(X.Fam[i,],X.Fam[j,])
   }
}

Relationship[1:5,1:5]

rel.pair <- Relationship[lower.tri(Relationship)]

colvec <- rep("blue",length(rel.pair))
colvec[rel.pair] <- "yellow"

plot(ibs.m,ibs.s,xlim=c(1,2),ylim=c(0,1),xlab="Mean",ylab="Standard deviation",
col=colvec)

legend(1.8,1,c("UN","PO"),col=c("blue","yellow"),pch=c(1,1))

#
# make (p0,p2) plot.
#

X[1:5,1:5]

stats <- function(x,y) {
  aux <- 2-abs(x-y) # number of shared alleles
  n0 <- sum(aux==0,na.rm=TRUE)
  n1 <- sum(aux==1,na.rm=TRUE)
  n2 <- sum(aux==2,na.rm=TRUE)
  n <- sum(!is.na(aux))
  p0 <- n0/n
  p1 <- n1/n
  p2 <- n2/n
  y <- c(p0,p1,p2)
  return(y)
}

stats(X[1,],X[2,])
sum(stats(X[1,],X[2,]))

Mp0 <- matrix(NA,nrow=n,ncol=n)
Mp1 <- matrix(NA,nrow=n,ncol=n)
Mp2 <- matrix(NA,nrow=n,ncol=n)

for(i in 1:n) {
   for(j in 1:n) {
      statsofapair <- stats(X[i,],X[j,])
      Mp0[i,j] <- statsofapair[1]
      Mp1[i,j] <- statsofapair[2]
      Mp2[i,j] <- statsofapair[3]
   }
}

p0vec <- Mp0[lower.tri(Mp0)]
p2vec <- Mp2[lower.tri(Mp2)]

plot(p0vec,p2vec,asp=1,col=colvec,
     xlab="% variants 0 shared alleles",
     ylab="% variants 2 shared alleles")

plot(jitter(p0vec,amount=.005),p2vec,asp=1,col=colvec,
xlab="% variants 0 shared alleles",
ylab="% variants 2 shared alleles")
legend("topright",c("UN","PO"),col=c("blue","yellow"),pch=c(1,1))

#
# IBD estimation with plink
#

getwd()
setwd("c:/plink")
getwd()

runstring <- "plink"
system(runstring)

runstring <- "plink --bfile YRI6 --genome --genome-full --out YRI6"
runstring
system(runstring)

Z <- read.table("YRI6.genome",header=TRUE)
colnames(Z)
head(Z)

table(rel.pair,Z$RT)

plot(Z$Z0,Z$Z1,asp=1,col=colvec,xlab=expression(k[0]),ylab=expression(k[1]))
legend("topright",c("UN","PO"),col=c("blue","yellow"),pch=c(1,1))

plot(jitter(Z$Z0,amount=0.05),Z$Z1,asp=1,col=colvec,xlab=expression(k[0]),ylab=expression(k[1]))
legend("topright",c("UN","PO"),col=c("blue","yellow"),pch=c(1,1))

