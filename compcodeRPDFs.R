# compcodeRPDFs.R
# R version 3.2.2 (2015-08-14)
# November 6, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating PDFs for the proportion of accurately identified genes from baySeq
# using simulated data from compcodeR. Files used in this script are the 
# output of file "compcodeRbaySeq.R".
#------------------------------------------------------------------------------
library(moments)

#setwd("/Users/mblai/Documents/GitHub/Homework-Repo")

# Read in input for the Negative Binomial (NB) prior and the Zero-Inflated Negative
# Binomial (ZINB) prior. 
NB <- read.csv(file = "compcodebaySeqOUT.csv")
ZINB <- read.csv(file = "compcodebaySeqZINBout.csv")

# Create summary statistics for each distribution.
# NB:
mean(NB$Correct.DE)
sd(NB$Correct.DE)

# Create a filtered PDF with varying filter widths for the NB prior.
plot(density(NB$Correct.DE, width = .01), ylab = 'PDF', xlab = 'x', main = '', type = 'n')
lines(density(NB$Correct.DE, width = .01), ylab = 'PDF', xlab = 'x', main = '')
lines(density(NB$Correct.DE, width = .03), ylab = 'PDF', xlab = 'x', main = '')


plot(density(ZINB$Correct.DE, width = .01), ylab = 'PDF', xlab = 'Correctly Identified Genes', main = '', type = "n")
lines(density(ZINB$Correct.DE, width = .01), ylab = 'PDF', xlab = 'Correctly Identified Genes', main = '')
lines(density(ZINB$Correct.DE, width = .03), ylab = 'PDF', xlab = 'Correctly Identified Genes', main = '')
lines(density(ZINB$Correct.DE, width = .05), ylab = 'PDF', xlab = 'Correctly Identified Genes', main = '')

xbar<- mean(ZINB$Correct.DE) 
xvar <- var(ZINB$Correct.DE)

a <- xbar*(((xbar*(1-xbar))/xvar)-1)
b <- (1-xbar)*(((xbar*(1-xbar))/xvar)-1)  

a/(a+b)

x <- seq(0, 1, by = .001)
lines(x, dbeta(x, shape1 = a, shape2 = b))
abline(v=xbar)

plot(density(ZINB$Correct.DE))

lines(density(ZINB$Correct.DE, width = .3), ylab = 'f', xlab = 'x', main = '')
plot(density(ZINB$Correct.DE))


kurtosis(ZINB$Correct.DE)

mean(ZINB$Correct.DE)
mean(NB$Correct.DE)

var(ZINB$Correct.DE)
var(NB$Correct.DE)

nbzinb <- c(NB$Correct.DE, ZINB$Correct.DE)
plot(density(nbzinb))



xbar<- mean(nbzinb) 
xvar <- var(nbzinb)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

xmode <- getmode(nbzinb)

b <- (1 - xbar)*(((xbar*(1 - xbar))/xvar)-1)
a <- (xmode*b - 2*xmode + 1)/(1 - xmode) 
#a <- xbar*(((xbar*(1-xbar))/xvar)-1)
  
x <- seq(0, 1, by = .001)
lines(x, dbeta(x, shape1 = a, shape2 = b))
abline(v=xbar)

ZINB200 <- read.csv(file.choose())
lines(density(zinb), col = 'blue')
zinb <- c(zinb, NB$Correct.DE)

