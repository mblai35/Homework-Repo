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
#-----------------------------------------------------------------------

# Read in the data. 
dat<-read.csv(file = "/Users/mblai/Documents/Classes/Ewers_readings/160303_gsAdataForBrent.csv")

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
plot(gs, A)

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

writeLines(modelString, con = "TEMPmodel.txt")

# Create vector of parameters. 
parameters<-c("beta0", "beta1", "a", "sigma")

m1<- jags.model(file="TEMPmodel.txt",
          data = m1.dat,
          n.chains = 3,
          n.adapt = 500)

update( m1 , n.iter=500 )
codaSamples = coda.samples( m1 , 
                            variable.names=c("beta0","beta1", "sigma","a","A") ,
                            n.iter=10000 )

parameterNames = varnames(codaSamples) # get all parameter names
#for ( parName in parameterNames ) {
#  diagMCMC( codaObject=codaSamples , parName=parName ) 
#}

mcmcMat = as.matrix(codaSamples)
chainLength = nrow(mcmcMat)
stepVec = floor(seq(1,chainLength,length=20)) # subset of chain
Nvar = ncol(gs)

png("Model_1_codaSamples.png")
plot(codaSamples)
dev.off()

install.packages("mcmcplots")
library(mcmcplots)
mcmcplot(codaSamples)
summary(codaSamples)

HPDinterval(codaSamples, prob = .95)

png("Model_1_caterplot.png")
caterplot(codaSamples)
dev.off()

# Find the xProbe and predicted y columns:
xPcols = grep( "A" , colnames(mcmcMat) , value=FALSE )
yPcols = grep( "a" , colnames(mcmcMat) , value=FALSE )
# Find the extreme predicted values for graph axis limits:
xLim = quantile( mcmcMat[,yPcols] , probs=c(0.005,0.995) )
# Make the plots of the posterior predicted values:
for ( i in 1:length(xPcols) ) {
  openGraph(width=4,height=3)
  plotPost( mcmcMat[,yPcols[i]] , xlab="y" , xlim=xLim , cenTend="mean" ,
            main=bquote( "Posterior Predicted y for x = "
                         * .(mcmcMat[1,xPcols[i]]) )  )
}

dic.samples(m1, n.iter=1000, type = "pD")
source("DBDA2E-utilities.R")

plot(gs, A, ylab="Net Assimilation", xlab="Stomatal Conductance")
########################### Model 2 ####################################
m1.dat<- list(gs=log10(gs), A=A)
m1.inits<-list(list("beta0"=1,"beta1"=-2, "tau"= 2)) 

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

writeLines(modelString, con = "TEMPmodel.txt")

parameters<-c("beta0", "beta1", "a", "sigma")

library(rjags)
library(coda)

m1<- jags.model(file="TEMPmodel.txt",
                data = m1.dat,
                n.chains = 3,
                n.adapt = 500)

update( m1 , n.iter=500 )
codaSamples = coda.samples( m1 , 
                            variable.names=c("beta0","beta1", "sigma","a","A") ,
                            n.iter=10000 )

parameterNames = varnames(codaSamples) # get all parameter names
#for ( parName in parameterNames ) {
#  diagMCMC( codaObject=codaSamples , parName=parName ) 
#}

mcmcMat = as.matrix(codaSamples)
chainLength = nrow(mcmcMat)
stepVec = floor(seq(1,chainLength,length=20)) # subset of chain
Nvar = ncol(gs)

dic.samples(m1, n.iter=1000, type = "pD")
########################### Model 3 ####################################

m1.dat<- list(gs=sqrt(gs), A=A)
m1.inits<-list(list("beta0"=1,"beta1"=-2, "tau"= 2)) 

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

writeLines(modelString, con = "TEMPmodel.txt")

parameters<-c("beta0", "beta1", "a", "sigma")

library(rjags)
library(coda)

m1<- jags.model(file="TEMPmodel.txt",
                data = m1.dat,
                n.chains = 3,
                n.adapt = 500)

update( m1 , n.iter=500 )
codaSamples = coda.samples( m1 , 
                            variable.names=c("beta0","beta1", "sigma","a","A") ,
                            n.iter=10000 )

parameterNames = varnames(codaSamples) # get all parameter names
#for ( parName in parameterNames ) {
#  diagMCMC( codaObject=codaSamples , parName=parName ) 
#}

mcmcMat = as.matrix(codaSamples)
chainLength = nrow(mcmcMat)
stepVec = floor(seq(1,chainLength,length=20)) # subset of chain
Nvar = ncol(gs)

png("Model_3_codaSamples.png",
    width=1700, height=1950,res=300)
plot(codaSamples[ ,c("beta0", "beta1")])
dev.off()

summary(codaSamples)
hpd<-HPDinterval(codaSamples, prob = .95)

png("Model_3_withCI.png",
    width=2000, height=2000,res=300)
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

