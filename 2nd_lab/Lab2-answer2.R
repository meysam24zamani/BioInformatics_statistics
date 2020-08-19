
##1. The file TSIChr22v4.raw contains genotype information of individuals from Tuscany in Italy,
##taken from the 1,000 Genomes project. The datafile contains all single nucleotide polymorphisms on 
##chromosome 22 for which complete information is available. Load this data into the R environment. 
##Use the fread instruction of the package data.table, which is more efficient for reading large datafiles. 
##This data is in (0,1,2) format, where 0 and 2 represent the homozygotes AA and BB, and 1 represents the heterozygote AB.
##The first six leading columns of the data matrix can be ignored, as they do not contain any genetic information.

library(data.table)
download.file("http://www-eio.upc.es/~jan/data/bsg/TSIChr22v4.raw","TSIChr22v4.raw")
data <- fread("TSIChr22v4.raw")



##2. How many individuals does the database contain, and how many variants? What percentage of the variants is monomorphic? 
##Remove all monomorphic SNPs from the database. How many variants remain in the database?

dim(data)

n <- nrow(data)
p <- ncol(data)
n

data <- data[,7:p]

monomorphic.AA <- which(apply(data,2,function(col) all(col==0)))
monomorphic.BB <- which(apply(data,2,function(col) all(col==2)))
monomorphics <- length(monomorphic.AA) + length(monomorphic.BB)
percentage.monomorphics <- (monomorphics/p)*100
X <- data[,-c(monomorphic.AA,monomorphic.BB)]
dim(X)

##3. Extract polymorphism rs587756191 T from the datamatrix, and determine its genotype counts.
##Apply a chi-square test for Hardy-Weinberg equilibrium, with and without continuity correction.
##Also try an exact test, and a permutation test. You can use function HWChisq, HWExact and HWPerm
##for this purpose. Do you think this variant is in equilibrium? Argue your answer.
genotypeOf <- genotype(X$rs587756191_T,sep="")
genotypeOfInfo <- summary(genotypeOf)
genotypeCounts <- genotypeOfInfo$genotype.freq[,1]
genotypeCounts
#The genotype counts are "genotypeCounts " 

#Xi-square test without continuity correction
chiNoCont <- HWChisq(X, cc= 0)
chiNoCont

#Xi-square test with continuity correction
chiCount <- HWChisq(X)
chiCount

#Exact test
exact <- HWExact(X)
exact

#Permutation test
permut <- HWPerm(X)
permut



##4. Determine the genotype counts for all these variants, and store them in a p x 3 matrix.**
##5. Apply a chi-square test without continuity correction for Hardy-Weinberg equilibrium to each SNP.
##You can use HWChisqStats for this purpose. How many SNPs are significant (use alpha = 0.05)
##Answer 4 & 5:
counts <- matrix(t(apply(X,2,function(col) c(sum(col==0),sum(col==1),sum(col==2)))), nrow=ncol(X), ncol=3)
colnames(counts) <- c("AA", "AB", "BB")
pvalues <- apply(counts, 1, function(row) HWChisq(row,cc = 0,verbose=FALSE)$pval)
passTest <- pvalues>= 0.05
sum(passTest)