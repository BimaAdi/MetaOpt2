# Shuffled Frog Leaping -----

source('./R/metaheuristic.FunctionCollection.R')

SFL <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
                numMemeplex=3, frogLeapingIteration=10){
  if(numMemeplex > numPopulation){
    stop("numMemeplex must less than or equal to numPopulation")
  }
  if(numMemeplex < 1){
    stop("numMemeplex must greater than 0")
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
  # initialize candidate solution
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineSFL(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution,
                       numMemeplex, frogLeapingIteration)
  return(bestPos)
}

engineSFL <- function(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution,
                      numMemeplex, frogLeapingIteration){
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  candidateSolutions <- data.frame(candidateSolution, fitness)
  
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # Sort candidate solutions based on fitness
    candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness),]
    
    # Split candidate solution into "numMemeplex" matrices
    memeplexId <- rep(1:numMemeplex,  length.out = numPopulation)
    candidateSolutions$memeplexId <- memeplexId
    
    # Memeplex step
    for(index in 1:numMemeplex){
      # Take memeplex with memeplexid = index
      memeplex <- candidateSolutions[candidateSolutions$memeplexId == index, ]
      for(i in 1:frogLeapingIteration){
        # Determine best and worst fitness in this memeplex
        memeplex <- memeplex[order(memeplex$fitness),]
        best <- as.matrix(head(memeplex[, 1:numVar], 1))
        worst <- as.matrix(tail(memeplex[, 1:numVar], 1))
        randomMatrix <- apply(matrix(rep(NA, numVar), nrow = 1, byrow = TRUE), c(1, 2), function(x){
          runif(1)
        })
        # Update worst using
        new <- worst + randomMatrix * (best - worst)
        # if worst have better fitness
        if(calcFitness(FUN, optimType, new) > calcFitness(FUN, optimType, worst)){
          new <- worst 
        }
        # Update memeplex
        memeplex[nrow(memeplex), 1:numVar] <- new
        memeplex$fitness <- calcFitness(FUN, optimType, as.matrix(memeplex[,1:numVar]))
      }
      # Shuffled step
      candidateSolutions[candidateSolutions$memeplexId == index, ] <- memeplex
    }
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, -1*optimType, as.matrix(candidateSolutions[, 1:numVar])))
}
