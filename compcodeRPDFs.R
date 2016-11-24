# compcodeRPDFs.R
# R version 3.2.2 (2015-08-14)
# November 6, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating PDFs for the proportion of accurately identified genes from baySeq
# using simulated data from compcodeR. Files used in this script are the 
# output of file "compcodeRbaySeq.R".
#------------------------------------------------------------------------------
library(moments)

# Read in input for the Negative Binomial (NB) prior and the Zero-Inflated Negative
# Binomial (ZINB) prior. 
NB <- read.csv(file.choose())
ZINB <- read.csv(file.choose())

# Create a filtered PDF with filter width of .035 for each prior.
plot(density(NB$Correct.DE, width = .035), ylab = 'PDF', xlab = 'Proportion Correct', 
     main = 'Correctly Identified Genes', type = 'n')
lines(density(NB$Correct.DE, width = .035))
lines(density(ZINB$Correct.DE, width = .035), lty = 2)
legend('topleft', c('NB', 'ZINB'), lty = c(1,2))

# Calculate statistics for each prior. 
# mean
mean(ZINB$Correct.DE)
mean(NB$Correct.DE)
# variance
var(ZINB$Correct.DE)
var(NB$Correct.DE)

# Merge NB and ZINB data. 
total <- rbind(NB, ZINB)

# Plot filtered PDF for combined NB and ZINB data.
plot(density(total$Correct.DE, width = .035), col = 'red')

# Find alpha and beta for Beta distribution from mean and variance.
xbar<- mean(total$Correct.DE) 
xvar <- var(total$Correct.DE)

a <- xbar*(((xbar*(1-xbar))/xvar)-1)
b <- (1-xbar)*(((xbar*(1-xbar))/xvar)-1)  

# Plot the Beta distribution with the appropriate alpha and beta.
x <- seq(0, 1, by = .001)
lines(x, dbeta(x, shape1 = a, shape2 = b))

# Find mode. 
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

xmode <- getmode(total$Correct.DE)

b <- (1 - xbar)*(((xbar*(1 - xbar))/xvar)-1)
a <- (xmode*b - 2*xmode + 1)/(1 - xmode) 
#a <- xbar*(((xbar*(1-xbar))/xvar)-1)

x <- seq(0, 1, by = .001)
lines(x, dbeta(x, shape1 = a, shape2 = b), lty =2)
abline(v=xbar)


