#1. The file rs394221.dat contains genotype information, for cases and controls, of polymorphism rs394221,
#which is presumably related to Alzheimer's disease.Load the data file into the R environment.
initialDataTogether <- read.delim("rs394221.dat", header = FALSE)

#2. (1p) What is the sample size? What is the number of cases and the number of controls? Construct the
#contingency table of genotype by case/control status.

summaryInTable <- do.call(cbind, lapply(initialDataTogether, summary))
Cases <- c(summaryInTable[3],summaryInTable[2],summaryInTable[1])
Controls <- c(summaryInTable[6],summaryInTable[5],summaryInTable[4])
tableX <- rbind(Cases,Controls)
colnames(tableX) <- c("MM","Mm","mm")

#The sample size is 1167
sum(tableX)

#The number of cases are 509 and the number of controls are 658
rowSums(tableX)

#The contingency is the following:
tableX

#3. (1p) Explore the data by plotting the percentage of cases as a function of the genotype, ordering the latter
#according to the number of M alleles. Which allele increases the risk of the disease?

totalGenotypeCounts <- colSums(tableX)
riskOfDisease <- Cases/totalGenotypeCounts

plot(c(0,1,2),c(riskOfDisease[3],riskOfDisease[2],riskOfDisease[1]),ylim=c(0,1),type="b", xlab="Genotype",ylab="Risk")

#From the plot we can see the allele that much increase the risk of desease seems to be M

#4. (2p) Test for equality of allele frequencies in cases and controls by doing an alleles test. Report the
#test statistic, its reference distribution, and the p-value of the test. Is there evidence for different allele
#frequencies?

tableY <- cbind(2*tableX[,1]+tableX[,2],2*tableX[,3]+tableX[,2])
colnames(tableY) <- c("M","m")

resultsChisqTest <- chisq.test(tableY,correct=FALSE)

#Test statistic and reference distribution:
resultsChisqTest
resultsChisqTest$expected

#The p-value of the test is 0.0002037, as it is very low and less than the significance level, we cannot accept the null
#hypotesis, concluding that there are evidences for differences of allele frequencies.

resultsChisqTest$p.value

#5. (2p) Which are the assumptions made by the alleles test? Perform and report any addtional tests you
#consider adequate to verify the assumptions. Do you think the assumptions of the alleles test are met?

#The test for equality of allele frequencies assumes independence.

#We will perform the fisher exact test to verify the assumptions.
fisher.test(tableY)
tableY

OR <- (tableY[1,1]*tableY[2,2])/(tableY[1,2]*tableY[2,1])

#We can check assumptions of the alleles tests are not met because we know OR = 1.364592 > 1 which indicates
#association.
OR

seLnOr <- sqrt(sum(1/tableY))
seLnOr
llLogOdds <- log(OR) - qnorm(0.975)*seLnOr
ulLogOdds <- log(OR) + qnorm(0.975)*seLnOr
llOdds <- exp(llLogOdds)
ulOdds <- exp(ulLogOdds)

#We also observed that the presence of m raises the odds of M, and the presence of M raises the odds of m.
llOdds
ulOdds


