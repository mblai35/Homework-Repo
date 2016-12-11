# MCMCproject.R
# R version 3.2.2 (2015-08-14)
# November 26, 2016. Mallory B. Lai.
# Code to accompany MCMC ppt presentation. 
#-----------------------------------------------------------------------
library(markovchain)
library(ggplot2)
#-----------------------------------------------------------------------

###### Metropolis Algorithm applied to States example:
set.seed(3)
# Number of steps for algorithm to make.
steps <- 100
# Initial position. 
current <- 3
# Initialize vector to store moves. 
movesMatrix <- matrix(rep(0, steps*2), ncol = 2)
# Name movesMatrix columns.
colnames(movesMatrix) <- c('TimeStep', 'State')

# Metropolis algorithm:
for (i in 1:steps)
{
  # Record current time step.
  movesMatrix[i, 1] <- i
  
  # Record current position.
  movesMatrix[i, 2] <- current
  
  # Simulate proposal to left or right with coin flip.
  proposal <- current + sample(c(-1, 1), size = 1)
  # Account for edges:
  if (proposal < 1) proposal <- 1
  if (proposal > 5) proposal <- 5
  
  # Probability of moving to proposed position.
  prob <- proposal/current
  current <- ifelse(runif(1) < prob, proposal, current)
  
}

#### Plot movement over time. 
# Convert matrix to data.frame for ggplot2.
movesMatrix <- data.frame(movesMatrix)
# Plot TimeStep as y-axis and State as x-axis and store as plot1.
(plot1 <- ggplot(movesMatrix, aes(State, TimeStep)) + geom_path())
# Plot bar chart of states.
(qplot1 <- qplot(State, data = movesMatrix, geom = "bar",
      xlab = "State", fill = factor(State)))

### Write plots to file for presentation: 
# Set dimensions for plot1 TIFF file.
tiff("MovesStatesEx1.tiff", width = 4, height = 8, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
plot1
# Reset plotting environment. 
dev.off()

# Set dimensions for qplot1 TIFF file.
tiff("BarChartStatesEx1.tiff", width = 4, height = 3, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
qplot1
# Reset plotting environment. 
dev.off()


### Repeat for another chain:
# Metropolis Algorithm:
set.seed(4)
# Number of steps for algorithm to make.
steps <- 100
# Initial position. 
current <- 3
# Initialize vector to store moves. 
movesMatrix$Chain2 <- c(rep(0, steps))

# Metropolis algorithm:
for (i in 1:steps)
{
  # Record current position.
  movesMatrix[i, 3] <- current
  
  # Simulate proposal to left or right with coin flip.
  proposal <- current + sample(c(-1, 1), size = 1)
  # Account for edges:
  if (proposal < 1) proposal <- 1
  if (proposal > 5) proposal <- 5
  
  # Probability of moving to proposed position.
  prob <- proposal/current
  current <- ifelse(runif(1) < prob, proposal, current)
  
}

# Plot bar chart of states.
(qplot2 <- qplot(Chain2, data = movesMatrix, geom = "bar",
                 xlab = "State", fill = factor(Chain2)))


# Plot movement over time for both chains.  
(plot2 <- plot1 + geom_path(data = movesMatrix, 
                           aes(Chain2, TimeStep, colour = "Chain2")))

### Write plots to file for presentation: 
# Set dimensions for plot2 TIFF file.
tiff("MovesStatesEx2.tiff", width = 4, height = 8, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
plot2
# Reset plotting environment. 
dev.off()

# Set dimensions for qplot1 TIFF file.
tiff("BarChartStatesEx2.tiff", width = 4, height = 3, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
qplot2
# Reset plotting environment. 
dev.off()

### Now up the number of steps:
set.seed(3)
# Number of steps for algorithm to make.
steps <- 1000
# Initial position. 
current <- 3
# Initialize vector to store moves. 
movesMatrix <- matrix(rep(0, steps*2), ncol = 2)
# Name movesMatrix columns.
colnames(movesMatrix) <- c('TimeStep', 'State')

# Metropolis algorithm:
for (i in 1:steps)
{
  # Record current time step.
  movesMatrix[i, 1] <- i
  
  # Record current position.
  movesMatrix[i, 2] <- current
  
  # Simulate proposal to left or right with coin flip.
  proposal <- current + sample(c(-1, 1), size = 1)
  # Account for edges:
  if (proposal < 1) proposal <- 1
  if (proposal > 5) proposal <- 5
  
  # Probability of moving to proposed position.
  prob <- proposal/current
  current <- ifelse(runif(1) < prob, proposal, current)
  
}

#### Plot movement over time. 
# Convert matrix to data.frame for ggplot2.
movesMatrix <- data.frame(movesMatrix)
# Plot TimeStep as y-axis and State as x-axis and store as plot1.
(plot3 <- ggplot(movesMatrix, aes(State, TimeStep)) + geom_path())
# Plot bar chart of states.
(qplot3 <- qplot(State, data = movesMatrix, geom = "bar",
                 xlab = "State", fill = factor(State)))

### Write plots to file for presentation: 
# Set dimensions for plot1 TIFF file.
#tiff("MovesStatesEx3.tiff", width = 4, height = 8, 
#     units = 'in', res = 300)
# Write plot to TIFF file. 
plot3
# Reset plotting environment. 
#dev.off()

# Set dimensions for qplot1 TIFF file.
#tiff("BarChartStatesEx3.tiff", width = 4, height = 3, 
#     units = 'in', res = 300)
# Write plot to TIFF file. 
qplot3
# Reset plotting environment. 
#dev.off()

## Second chain for 1000 steps
set.seed(4)
# Number of steps for algorithm to make.
steps <- 1000
# Initial position. 
current <- 3
# Initialize vector to store moves. 
movesMatrix$Chain2 <- c(rep(0, steps))

# Metropolis algorithm:
for (i in 1:steps)
{
  # Record current time step.
  #movesMatrix[i, 1] <- i
  
  # Record current position.
  movesMatrix[i, 3] <- current
  
  # Simulate proposal to left or right with coin flip.
  proposal <- current + sample(c(-1, 1), size = 1)
  # Account for edges:
  if (proposal < 1) proposal <- 1
  if (proposal > 5) proposal <- 5
  
  # Probability of moving to proposed position.
  prob <- proposal/current
  current <- ifelse(runif(1) < prob, proposal, current)
  
}

# Plot bar chart of states.
(qplot3.b <- qplot(Chain2, data = movesMatrix, geom = "bar",
                 xlab = "State", fill = factor(Chain2)))

# Plot movement over time for both chains.  
(plot3.b <- plot3 + geom_path(data = movesMatrix, 
                            aes(Chain2, TimeStep, colour = "Chain2")))

### Write plots to file for presentation: 
# Set dimensions for plot2 TIFF file.
tiff("MovesStatesEx3.b.tiff", width = 4, height = 8, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
plot3.b
# Reset plotting environment. 
dev.off()

# Set dimensions for qplot1 TIFF file.
tiff("BarChartStatesEx3.b.tiff", width = 4, height = 3, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
qplot3.b
# Reset plotting environment. 
dev.off()


### Now up the number of steps AGAIN:
set.seed(3)
# Number of steps for algorithm to make.
steps <- 10000
# Initial position. 
current <- 3
# Initialize vector to store moves. 
movesMatrix <- matrix(rep(0, steps*2), ncol = 2)
# Name movesMatrix columns.
colnames(movesMatrix) <- c('TimeStep', 'State')

# Metropolis algorithm:
for (i in 1:steps)
{
  # Record current time step.
  movesMatrix[i, 1] <- i
  
  # Record current position.
  movesMatrix[i, 2] <- current
  
  # Simulate proposal to left or right with coin flip.
  proposal <- current + sample(c(-1, 1), size = 1)
  # Account for edges:
  if (proposal < 1) proposal <- 1
  if (proposal > 5) proposal <- 5
  
  # Probability of moving to proposed position.
  prob <- proposal/current
  current <- ifelse(runif(1) < prob, proposal, current)
  
}

#### Plot movement over time. 
# Convert matrix to data.frame for ggplot2.
movesMatrix <- data.frame(movesMatrix)
# Plot TimeStep as y-axis and State as x-axis and store as plot1.
(plot4 <- ggplot(movesMatrix, aes(State, TimeStep)) + geom_path())
# Plot bar chart of states.
(qplot4 <- qplot(State, data = movesMatrix, geom = "bar",
                 xlab = "State", fill = factor(State)))

### Write plots to file for presentation: 
# Set dimensions for plot1 TIFF file.
tiff("MovesStatesEx4.tiff", width = 4, height = 8, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
plot4
# Reset plotting environment. 
dev.off()

# Set dimensions for qplot1 TIFF file.
tiff("BarChartStatesEx4.tiff", width = 4, height = 3, 
     units = 'in', res = 300)
# Write plot to TIFF file. 
qplot4
# Reset plotting environment. 
dev.off()





# Create transition matrix
tMatrix <- matrix(c(.5, .5, 0, 0, 0,
                    .25, .25, .5, 0, 0,
                    0, 1/3, 1/6, .5, 0, 
                    0, 0, .375, .125, .5,
                    0, 0, 0, .4, .6), 
                  ncol = 5, byrow = T)

# Create a discrete time Markov Chain using the package "markovchain".
MC1 <- as(tMatrix, "markovchain")

# Set seed to prevent graph from changing. 
set.seed(2)
# Graph MC1. 
plot(MC1, vertex.size = 20)

# Examine summary of MC1
summary(MC1)

# Find the equilibrium distribution.
steadyStates(MC1)[1, ]

# Create matrix with states and steady state solutions for plotting.
# Multiply steady state solution by 15 to restore original values. 
SSmatrix <- data.frame(MC1@states, steadyStates(MC1)[1, ]*15)
colnames(SSmatrix) <- c('state', 'ss')
# Bar chart of population for each state. 
qplot(state, data = SSmatrix, geom = "bar", weight = ss, 
      ylab = "Population", xlab = "State", fill = state)
#ggplot(smat, aes(state, ss)) + geom_point(aes(size = ss))

