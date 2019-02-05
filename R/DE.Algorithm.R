# Differential Evolution (DE)

source('./R/metaheuristic.FunctionCollection.R')

FUN <- sphere
optimType <- "MIN"
numPopulation <- 10
numVar <- 2
maxIter <- 500
rangeVar <- matrix(c(-10, 10), nrow = 2)
scalingVector <- 0.1
crossOverRate <- 0.5

# calculate the dimension of problem if not specified by user
dimension <- ncol(rangeVar)

# parsing rangeVar to lowerBound and upperBound
lowerBound <- rangeVar[1,]
upperBound <- rangeVar[2,]

# if user define the same upper bound and lower bound for each dimension
if(dimension==1){
  dimension <- numVar
}

## convert optimType to numerical form
## 1 for minimization and -1 for maximization
if(optimType == "MAX") optimType <- -1 else optimType <- 1

# if user only define one lb and ub, then repeat it until the dimension
if(length(lowerBound)==1 & length(upperBound)==1){
  lowerBound <- rep(lowerBound, dimension)
  upperBound <- rep(upperBound, dimension)
}

# generate candidate solution
candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
fitness <- calcFitness(FUN, optimType, candidateSolution)
candidateSolutions <- data.frame(candidateSolution, fitness)

newCandidateSolutions <- candidateSolutions
for(row in 1:nrow(candidateSolutions)){
  # create mutant vector
  index <- c(1:numPopulation, 1, 2)#index to prevent overindexing
  mutant1 <- data.matrix(candidateSolutions[index[row], 1:numVar])
  mutant2 <- data.matrix(candidateSolutions[index[row + 1], 1:numVar])
  mutant3 <- data.matrix(candidateSolutions[index[row + 2], 1:numVar])
  mutantVector <- mutant1 + scalingVector*(mutant2 - mutant3)
  
  # crossover step
  crossoverVector <- c()
  for(variable in 1:length(mutantVector)){
    if(runif(1) <= crossOverRate){
      crossoverVector[variable] <- mutantVector[variable]
    }else{
      crossoverVector[variable] <- candidateSolutions[row, variable]
    }
  }
  
  # selection step
  crossoverVectorFitness <- calcFitness(FUN, optimType, matrix(crossoverVector, ncol = numVar, byrow = TRUE))
  previousCandidateSolutionFitness <- candidateSolutions[row, "fitness"]
  if(crossoverVectorFitness > previousCandidateSolutionFitness){
    candidateSolution <- data.matrix(candidateSolutions[row, 1:numVar])
    fitness <- previousCandidateSolutionFitness
    newCandidateSolutions[row, ] <- c(candidateSolution, fitness)
  }else{
    candidateSolution <- matrix(crossoverVector, ncol = numVar, byrow = TRUE)
    fitness <- crossoverVectorFitness
    newCandidateSolutions[row, ] <- c(candidateSolution, fitness)
  }
}
candidateSolutions <- newCandidateSolutions
calcBest(FUN, -1*optimType, as.matrix(candidateSolutions[, 1:numVar]))
a <- function(a=0.2){
  if(a < 0 || a > 1){
    #stop("a must between 0 and 1")
    warning("you fool")
  }
  return(a *10)
}
a <- matrix(rep(2, 50), ncol = 5, byrow = TRUE)
t(apply(a, 1, function(x){
  print(x)
  c(1, 2, 3, 4, 5)
}))