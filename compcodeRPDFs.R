# compcodeRPDFs.R
# R version 3.2.2 (2015-08-14)
# November 6, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Creating PDFs for the proportion of accurately identified genes 
# from baySeq using simulated data from compcodeR. Files used in 
# this script are the output of file "compcodeRbaySeq.R".
#-----------------------------------------------------------------------
library(knitr)

# Read in input for the Negative Binomial (NB) prior and the 
# Zero-Inflated Negative Binomial (ZINB) prior. 
NB <- read.csv(file.choose())
ZINB <- read.csv(file.choose())

#### Create a filtered PDF for DR with filter width of .035 for each 
#### prior.
# Initialize plot.
plot(density(ZINB$Correct.DE, width = .035), 
     ylab = 'PDF', xlab = 'Detection Rate',
     main = 'Proportion of True DE Genes Detected', 
     xlim = c(0.1, .35), type = 'n') 

# Add NB prior filtered PDF.
lines(density(NB$Correct.DE, width = .035)) 

# Add ZINB prior filtered PDF.
lines(density(ZINB$Correct.DE, width = .035), lty = 5) 

# Add legend.
legend('topleft', c('NB', 'ZINB'), lty = c(1,5)) 

#### Calculate DR statistics for each prior and put in a table.
# Create a matrix holding the mean and variance for the NB and 
# ZINB priors.
statistics <- matrix(c(mean(ZINB$Correct.DE), 
                       mean(NB$Correct.DE), var(ZINB$Correct.DE), 
                       var(NB$Correct.DE)), nrow = 2)

# Name the rows 'NB' and 'ZINB'.
rownames(statistics) <- c('NB', 'ZINB')

# Name the columns 'Mean' and 'Variance'. 
colnames(statistics) <- c('Mean', 'Variance')

# Create a table holding the priors' stats.
kable(statistics)

#### Fit Beta PDF to combined filtered PDF
# Merge NB and ZINB data. 
total <- rbind(NB, ZINB)

# Find alpha and beta for Beta distribution from mean and variance.
# Calculate mean.
xbar<- mean(total$Correct.DE) 

# Calculate variance.
xvar <- var(total$Correct.DE)

# Calculate alpha.
a <- xbar*(((xbar*(1-xbar))/xvar)-1)

# Calculate beta. 
b <- (1-xbar)*(((xbar*(1-xbar))/xvar)-1)  

#### Plot filtered PDF for combined NB and ZINB data with 
#### Beta distribution.
# Initialize a sequence of values between 0 and 1 to plot.
x <- seq(0, 1, by = .001)

# Initalize plot. 
plot(x, dbeta(x, shape1 = round(a), shape2 = round(b)), 
     xlab = 'Proportion', ylab = 'PDF', 
     main = 'True DE Genes Detected', type = 'n')

# Plot filtered PDF for combined NB and ZINB prior.
lines(density(total$Correct.DE, width = .035))

# Plot Beta PDF.
lines(x, dbeta(x, shape1 = round(a), shape2 = round(b)), lty = 2)

# Add legend.
legend('topright', c('Filtered PDF', 'Beta(46, 141)'), lty = c(1,2))

#### Find alpha and beta for Beta distribution from mean and 
#### variance for TPR.
# Calculate mean.
xbar<- mean(total$Correct.TopCounts) 

# Calculate variance. 
xvar <- var(total$Correct.TopCounts)

# Calculate alpha.
a <- xbar*(((xbar*(1-xbar))/xvar)-1)

# Calculate beta.
b <- (1-xbar)*(((xbar*(1-xbar))/xvar)-1)  

#### Plot the Beta distribution with the appropriate alpha and beta.
# Initialize a sequence of values between 0 and 1 to plot.
x <- seq(0, 1, by = .001)

# Plot Beta PDF for TPR.
plot(x, dbeta(x, shape1 = round(a), shape2 = round(b)), 
     lty = 2, xlim = c(.9, 1),
     ylim = c(0, 60), type = 'l', 
     ylab = "PDF", xlab = "Proportion", main ="True Positive Rate")

# Add Filtered PDF from combined NB and ZINB prior.
lines(density(total$Correct.TopCounts))

# Add legend.
legend('topleft', c('Filtered PDF', 'Beta(227, 5)'), lty = c(1,2)) 



