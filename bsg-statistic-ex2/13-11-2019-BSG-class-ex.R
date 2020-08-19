#Class exercise:

install.packages("HardyWeinberg")
library(HardyWeinberg)

x <- c(CC=23,CT=48,TT=29)
x
results <- HWChisq(x,cc=0)

results.cc <- HWChisq(x)

results.exact <- HWExact(x)


# Next question 

y <- c(CC=0,CT=7,TT=93)
y

resultsy <- HWChisq(y,cc=0)

resultsy.cc <- HWChisq(y)

resultsy.exact <- HWExact(y)

z <- rbind(x,y)
z

HWTernaryPlot(z)

HWExact(x)
?HWPerm
HWPerm(x,nperm = 50000)



HWPerm(y)

?HWData
set.seed(123)
z <- HWData(nm=100)
z

HWTernaryPlot(z)

?HWChisq

out <- HWChisqStats(z)
attributes(out)
out

cbind(z,out)

sum(out<0.05,na.rm=TRUE)


?HWExact

pval.exa <- HWExactStats(z)
sum(pval.exa < 0.05)


z <- HWData(nm=10000)
chi <- HWChisqStats(z)
chi
hist(chi)
max(chi,ns.rm=TRUE)
hist(chi,breaks = 50)
