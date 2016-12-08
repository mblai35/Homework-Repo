# ecophysMCMC.R
# R version 3.2.2 (2015-08-14)
# November 26, 2016. Mallory B. Lai.
# Gibbs sampling applied to ecophysiological data from Brassica rapa. 
#-----------------------------------------------------------------------
#setwd("/Users/mblai/Downloads/DBDA2Eprograms")
#source("DBDA2E-utilities.R")
library(rjags)
library(runjags)
library(coda)
library(ggplot2)
#-----------------------------------------------------------------------

# Read in the data. 
dat <- read.csv(file = "/Users/mblai/Documents/Classes/Ewers_readings/160303_gsAdataForBrent.csv")

# Store stomatal conductance for all non-RIL genotypes.
gs <- c(dat$gs_.mmol.m.2.s.1.[dat$genotype == "bro"], 
             dat$gs_.mmol.m.2.s.1.[dat$genotype == "oil"],
             dat$gs_.mmol.m.2.s.1.[dat$genotype == "cab"],
             dat$gs_.mmol.m.2.s.1.[dat$genotype == "tur"])

# Store net assimilation for all non-RIL genotypes. 
A <- c(dat$A_.umol.m.2.s.1.[dat$genotype == "bro"], 
            dat$A_.umol.m.2.s.1.[dat$genotype == "oil"],
            dat$A_.umol.m.2.s.1.[dat$genotype == "cab"],
            dat$A_.umol.m.2.s.1.[dat$genotype == "tur"])

# Look at relationship between stomatal conductance and net assimilation
# with a basic plot. 
# Combine data into dataframe for ggplot2.
combDat <- data.frame(A, gs)

# Initialize a png file. 
png("A_vs_gs.png",
    width=2000, height=2000,res=300)
# Plot data.
ggplot(combDat, aes(gs, A)) + geom_point() + 
  labs(x = "Stomatal Conductance", y = "Net Assimilation",
  title = "Brassica crop accessions")
# Reset plotting environment. 
dev.off()

# It is evident from the plot that some sort of linear model will be
# an appropriate model to describe the data. Three different models
# will be compared: linear, log-transformed, and sqrt-transformed 
# models. 

########################### Model 1 ####################################

# Create a list to store data for a linear model. 
m1.dat <- list(gs=gs, A=A)
# Give the model initial values for the markov chain to start at.
m1.inits <- list(list("beta0" = 1,"beta1" = -2, "tau" = 2)) 

# Store model for rjags. 
modelString =
"model{
for(i in 1:length(gs))
{
  A[i] ~ dnorm(a[i], tau)
  a[i] <- beta0 + beta1 * gs[i]
}
beta0 ~ dnorm(0, 1*10^(-6))
beta1 ~ dnorm(0, 1*10^(-6))
tau ~ dgamma(.001,.001)
sigma <- 1/sqrt(tau)
}"

# Save model to a text file for rjags.
writeLines(modelString, con = "TEMPmodel.txt")

# Create vector of parameters. 
parameters <- c("beta0", "beta1", "a", "sigma")

# Compile model.
m1 <- jags.model(file = "TEMPmodel.txt",
          data = m1.dat,
          n.chains = 3,
          n.adapt = 500)

# Run model and store as mcmc.list object.
codaSamples = coda.samples(m1, 
              variable.names = c("beta0", "beta1", "sigma", "a", "A"),
              n.iter = 10000 )


#png("Model_1_codaSamples.png")
plot(codaSamples[ ,c("beta0", "beta1", "sigma")], 
     col = c('black', 'gold2', 'dodgerblue3'))
#dev.off()

# Find credible interval of beta0, beta1, and sigma.
HPDinterval(codaSamples[ ,c("beta0", "beta1", "sigma")], prob = .95)

# Calculate DIC of model 1 for model comparison. 
dic.samples(m1, n.iter=1000, type = "pD") # pD = penalized deviance

########################### Model 2 ####################################

# Create a list to store data for a log transformed linear model. 
m1.dat<- list(gs=log10(gs), A=A)

# Give the model initial values for the markov chain to start at.
m1.inits<-list(list("beta0"=1,"beta1"=-2, "tau"= 2)) 

# Store model for rjags.
modelString=
  "model{
for(i in 1:length(gs))
{
  A[i]~dnorm(a[i], tau)
  a[i]<-beta0+beta1*(gs[i])
}
beta0 ~dnorm(0, 1*10^(-6))
beta1 ~dnorm(0, 1*10^(-6))
tau ~dgamma(.001,.001)
sigma<- 1/sqrt(tau)
}"

# Write model to temporary text file for rjags to read.
writeLines(modelString, con = "TEMPmodel.txt")

# Compile model. 
m1<- jags.model(file="TEMPmodel.txt",
                data = m1.dat,
                n.chains = 3,
                n.adapt = 500)

# Run model and store as mcmc.list object. 
codaSamples = coda.samples( m1 , 
                            variable.names=c("beta0","beta1", "sigma","a","A") ,
                            n.iter=10000 )

# Plot trace plot. 
plot(codaSamples[ ,c("beta0", "beta1", "sigma")], 
     col = c('black', 'gold2', 'dodgerblue3'))

# Calculate DIC of model 2 for model comparison.
dic.samples(m1, n.iter=1000, type = "pD")

########################### Model 3 ####################################

# Create a list to store data for a log transformed linear model. 
m1.dat <- list(gs=sqrt(gs), A=A)

# Give the model initial values for the markov chain to start at.
m1.inits <- list(list("beta0" = 1,"beta1" = 50, "tau" = 1)) 

# Store model for rjags.
modelString=
  "model{
for(i in 1:length(gs))
{
  A[i]~dnorm(a[i], tau)
  a[i]<-beta0+beta1*gs[i]
}
beta0 ~dnorm(0, 1*10^(-6))
beta1 ~dnorm(0, 1*10^(-6))
tau ~dgamma(.001,.001)
sigma<- 1/sqrt(tau)
}"

# Write model to temporary text file for rjags to read.
writeLines(modelString, con = "TEMPmodel.txt")

# Compile the model. 
m1<- jags.model(file="TEMPmodel.txt",
                data = m1.dat,
                n.chains = 3,
                n.adapt = 500)

# Burn-in:
update(m1, n.iter = 1000)

# Run model and store as mcmc.list object. 
codaSamples = coda.samples( m1 , 
                n.adapt = 1000, 
                variable.names = c("beta0", "beta1", "sigma", "a", "A"),
                n.iter=10000 )

png("Model_3_codaSamplesBurnIn.png",
    width=1700, height=1950,res=300)
plot(codaSamples[ ,c("beta0", "beta1")], 
     col = c('black', 'gold2', 'dodgerblue3'))
dev.off()

# Calculate DIC of model 3 for model comparison.
dic.samples(m1, n.iter=1000, type = "pD")

# Find credible interval of beta0, beta1, and sigma.
HPDinterval(codaSamples[ ,c("beta0", "beta1", "sigma")], prob = .95)

# Since model 3 has the lowest DIC (although not considerably lower
# than model 2), we'll 
plot(gs, A, ylab="Net Assimilation", 
     xlab="Stomatal Conductance",
     main = "Square-Root Transformed 
    Bayesian Model with 95% CI")
curve(-5.5255+42.7834561*sqrt(x), add=T, lty=2, col="red")
curve(2.78+59.505*sqrt(x), add=T, lty=2, col="red")
curve(-1.39+ 51.1512*sqrt(x), add=T)
dev.off()

png("Model_3_predictionCI.png",
    width=2000, height=2000,res=300)
plot(gs, A, ylab="Net Assimilation", 
     xlab="Stomatal Conductance",
     main = "RIL Genotype with 95% Prediction CI",
     type="n")
curve(-5.5255+42.7834561*sqrt(x), add=T, lty=2, col="red")
curve(2.78+59.505*sqrt(x), add=T, lty=2, col="red")
curve(-1.39+ 51.1512*sqrt(x), add=T)
points(c(dat$gs_.mmol.m.2.s.1.[dat$genotype=="ril301"],
dat$gs_.mmol.m.2.s.1.[dat$genotype=="ril46"]), c(dat$A_.umol.m.2.s.1.[dat$genotype=="ril301"],
dat$A_.umol.m.2.s.1.[dat$genotype=="ril46"]),
col="black", pch=20)
dev.off()


# Initialize a png file. 
png("A_vs_gs.png",
    width=2000, height=2000,res=300)
# Plot data.
p1 <- ggplot(combDat, aes(gs, A)) + geom_point() + 
  labs(x = "Stomatal Conductance", y = "Net Assimilation",
       title = "Brassica crop accessions")
# Reset plotting environment. 
dev.off()

ggplot(combDat, aes(gs, A)) + geom_line() 

curve(-5.5255+42.7834561*sqrt(x), add=T, lty=2, col="red")
curve(2.78+59.505*sqrt(x), add=T, lty=2, col="red")
curve(-1.39+ 51.1512*sqrt(x), add=T)



