library(markovchain)
library(ggplot2)

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

# Metropolis Algorithm:
set.seed(3)
# Number of steps for algorithm to make.
steps <- 100
# Initial position. 
current <- 3
# Initialize vector to store moves. 
movesMatrix <- matrix(rep(0, steps*2), ncol = 2)
# Name movesMatrix columns.
colnames(movesMatrix) <- c('TimeStep', 'State')

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
plot1 <- ggplot(movesMatrix, aes(TimeStep, State)) + geom_path()
# Plot bar chart of states.
ggplot(movesMatrix, aes(State)) + geom_bar(aes(fill = State))

# Repeat for another chain:
# Metropolis Algorithm:
set.seed(4)
# Number of steps for algorithm to make.
steps <- 100
# Initial position. 
current <- 3
# Initialize vector to store moves. 
movesMatrix$State2 <- c(rep(0, steps))

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
ggplot(movesMatrix, aes(State2)) + geom_bar(aes(fill = State2))

# Plot movement over time for both chains.  
plot2 <- plot1 + geom_path(data = movesMatrix, aes(TimeStep,State2, 
                                                   colour = "State2"))




