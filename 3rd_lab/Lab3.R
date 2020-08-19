#install.packages("genetics")
#install.packages("HardyWeinberg")
#install.packages("snpStats")
#install.packages("LDheatmap")
library(genetics)
library(HardyWeinberg)
library(LDheatmap)
#2.Load the FOXP2.dat file into the R environment. How many individuals and how many SNPs
#are there in the database? What percentage of the data is missing?
data <- read.delim("FOXP2.dat", header = TRUE, row.names = 1, sep=" ")

removeSlashes <- function(x) {
  gsub("/", "", x, fixed=TRUE)
}
data <- as.data.frame(apply(data, c(1,2), removeSlashes))

nrows <- nrow(data)
ncols <- ncol(data)
# There are 104 individuals and 543 SNPs.

percentageMissing <- 100*sum(is.na(data))/(nrows*ncols)
# There are not missing data in this dataset.

#3. Determine the genotype counts for each SNP, and depict all SNPs simultaeneously in a ternary
#plot, and comment on your result. For how many variants do you reject Hardy-Weinberg equilibrium
#using an ordinary chi-square test without continuity correction? (hint: you can read the .bim in R
#in order to determine the alleles of each SNP, and use function MakeCounts from the HardyWeinberg
#package to create a matrix of genotype counts).

bimData <- read.delim("FOXP2.bim", header = FALSE, sep="\t")
alleles <- c(do.call(paste,c(bimData[c("V5", "V6")],sep="/")))
counts <- MakeCounts(data, alleles)

counts <- counts[,1:3]
results.chi <- apply(counts,1,HWChisq, verbose=FALSE, cc=0)

#Reject the ones with p-value > 0.05

filterReject <- function(x) {
  if (x$pval > 0.05) {
    return(x)
  }
}
filtered <- lapply(results.chi,filterReject)

filterNull <- function(x) {
  vector = c()
  for (i in 1:length(x)) {
    if (is.null(x[[i]])) {
      vector <- c(vector, i)  
    }
  }
  return(vector)
}
nullsIndex <- filterNull(filtered)
length(nullsIndex)
# We are going to reject 33 variants

#Ternary plot????????

#4. Using the function LD from the genetics package, compute the LD statistic D for the SNPs
#rs34684677 and rs2894715 of the database. Is there significant association between the alleles of
#these two SNPs?

datars34684677 <- data$rs34684677
datars2894715 <- data$rs2894715
gtrs34684677 <- genotype(datars34684677,sep="")
gtrs2894715 <- genotype(datars2894715,sep="")
ldQ4 <- LD(gtrs34684677,gtrs2894715)
ldQ4$D
ldQ4$`D'`
ldQ4$`P-value`

# The value of the LD statistic D is 0.9986536, which is very high, this means the imbalance between
#the alleles is big, which means the association between them are quite high. The high value of statistic
#D' also indicates the strong relation. Lastly, a P-value smaller than 0.05 indicates us the results are
#reliable.

#5. Also compute the LD statistic D for the SNPs rs34684677 and rs998302 of the database. Is
#there significant association between these two SNPs? Is there any reason why rs998302 could have
#stronger or weaker correlation than rs2894715?

datars998302 <- data$rs998302
gtrs998302 <- genotype(datars998302,sep="")
ldQ5 <- LD(gtrs34684677,gtrs998302)
ldQ5$D

# The value of the LD statistic D is 0.007208888, which is very low, this means the balance between
#the alleles is big, which means the association between them are small because they tends to be inherited
#randomly. The low value of statistic D' also indicates the weak relation. Lastly, a P-value smaller than 0.05
#indicates us the results are reliable.
# The reason for which rs2894715 could have stronger correlation than rs998302 can be that rs2894715
#is closer to rs34684677 in the chromosome.

#6. Given your previous estimate of D for SNPs rs34684677 and rs2894715, infer the haplotype
#frequencies. Which haplotype is the most common?
alleleFreq1 <- summary(gtrs34684677)
allele.freq[,1]
summary(gtrs2894715)

