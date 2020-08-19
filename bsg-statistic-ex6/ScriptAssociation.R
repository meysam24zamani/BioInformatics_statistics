#
# Script association analysis
#

Cases    <- c(112,278,150)
Controls <- c(206,348,150)

X <- rbind(Cases,Controls)
rownames(X) <- c("Cases","Controls")
colnames(X) <- c("MM","Mm","mm")

X
rowSums(X)
n <- sum(X)
n

Y <- cbind(2*X[,1]+X[,2],2*X[,3]+X[,2])
colnames(Y) <- c("M","m")
Y

sum(Y)
rowSums(Y)

#
# Alleles test
#

results <- chisq.test(Y,correct=FALSE)
results$expected

fisher.test(Y)

Y

or <- (Y[1,1]*Y[2,2])/(Y[1,2]*Y[2,1])
or

se.lor <- sqrt(sum(1/Y))
se.lor

ll.logodds <- log(or) - qnorm(0.975)*se.lor
ul.logodds <- log(or) + qnorm(0.975)*se.lor

ll.odds <- exp(ll.logodds)
ul.odds <- exp(ul.logodds)

ll.odds
ul.odds

#
# Armitage trend test
#

X

Cases
Controls
cas <- rep(c(0,1,2),Cases)
con <- rep(c(0,1,2),Controls)

x <- c(rep(1,sum(Cases)),
       rep(0,sum(Controls)))

y <- c(cas,con)

length(x)
length(y)

r <- cor(x,y)
r
n
A <- n*(r^2)
A

pvalue <- pchisq(A,df=1,lower.tail=FALSE)
pvalue

#
# Plot of risk as function of genotype
#

X

total.genotype.counts <- colSums(X)
total.genotype.counts

risk <- Cases/total.genotype.counts
risk

plot(c(0,1,2),risk,ylim=c(0,1),type="b",
xlab="Genotype",ylab="Risk")

lm1 <- lm(x~y)
summary(lm1)
abline(coefficients(lm1),lty="dotted",col="red")

#
# Co-dominant model
#

chisq.test(X)
fisher.test(X)

#
# Dominant model
#

X
Z <- cbind(X[,1],X[,2]+X[,3])
colnames(Z) <- c("MM","not MM")

chisq.test(Z,correct=FALSE) 

#
# Recessive model
#

W <- cbind(X[,1]+X[,2],X[,3])
colnames(W) <- c("not mm","mm")
chisq.test(W)

#
# Logistic regression
#

y
x

newy <- x
newx <- y
newx

#
# Genotype treated as quantitative
#

out.lm <- glm(newy~newx, 
          family = binomial(link = "logit"))


summary(out.lm)
b <- coefficients(out.lm)
b

OR <- exp(b[2])
OR

V <- vcov(out.lm)
se <- sqrt(diag(V))
se

ll <- b[2]-qnorm(0.975)*se[2]
ul <- b[2]+qnorm(0.975)*se[2]

ll.or <- exp(ll)
ul.or <- exp(ul)

OR
ll.or
ul.or

#
# Genotype as qualitative
#

newx

x.cat <- rep(NA,length(newx))
x.cat[newx==0] <- "MM"
x.cat[newx==1] <- "Mm"
x.cat[newx==2] <- "mm"

table(x.cat)

X

x.cat <- factor(x.cat)
x.cat

out1.lm <- glm(newy~x.cat, 
          family = binomial(link = "logit"))


summary(out1.lm)
b <- coefficients(out1.lm)
b

ORs <- exp(b)
ORs















