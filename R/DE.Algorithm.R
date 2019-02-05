# Differential Evolution (DE)

source('./R/metaheuristic.FunctionCollection.R')

DE <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
               scalingVector=0.1, crossOverRate=0.5){
  # check parameter scalingVector
  if(scalingVector < 0 || scalingVector > 1){
    stop("parameter scalingVector must between 0 and 1")
  }
  
  # check parameter crossOverRate
  if(crossOverRate < 0 || crossOverRate > 1){
    stop("parameter crossOverRate must between 0 and 1")
  }
  
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
  bestPos <- engineDE(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution,
                      scalingVector, crossOverRate)
  
  return(bestPos)
}

engineDE <- function(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution,
                     scalingVector, crossOverRate){
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  candidateSolutions <- data.frame(candidateSolution, fitness)
  
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
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
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, -1*optimType, as.matrix(candidateSolutions[, 1:numVar])))
}