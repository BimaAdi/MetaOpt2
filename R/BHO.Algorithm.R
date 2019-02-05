# Black Hole-based Optimization (BHO)

source('./R/metaheuristic.FunctionCollection.R')

BHO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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
  bestPos <- engineBHO(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution)
  
  return(bestPos)
}

engineBHO <- function(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution){
  # generate candidate solutions
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  candidateSolutions <- data.frame(candidateSolution, fitness)
  
  for(t in 1:maxIter){
    # select best candidate solution as blackhole and other as star
    bestCandidateSolution <- order(candidateSolutions$fitness)[1]
    blackhole <- candidateSolutions[bestCandidateSolution, ]
    star <- candidateSolutions[-1*bestCandidateSolution, ]
    
    # change location each star candidate solution (cs)
    csStar <- data.matrix(star[,1:numVar])
    randomMatrix <- matrix(runif(numVar * nrow(star)), ncol = numVar)
    csBlackhole <- matrix(rep(as.numeric(blackhole[1, 1:numVar]), nrow(star)), ncol = numVar, byrow = TRUE)
    star[,1:numVar] <- csStar + randomMatrix * (csBlackhole - csStar)
    star$fitness <- calcFitness(FUN, optimType, data.matrix(star[,1:numVar]))
    
    # if a star reaches better fitness exchange it with blackhole
    bestStar <- star[order(star$fitness)[1], ]
    if(bestStar$fitness < blackhole$fitness){
      temp <- blackhole
      blackhole <- bestStar
      star[order(star$fitness)[1], ] <- temp
    }
    
    # if a star cross event horizon generate new candidate solution for that star
    eventHorizon <- blackhole$fitness/sum(c(blackhole$fitness, star$fitness))
    isCrossEventHorizon <- abs(star$fitness - blackhole$fitness) < eventHorizon
    numStarThatCrossEventHorizon <- length(which(isCrossEventHorizon == TRUE))
    if(numStarThatCrossEventHorizon > 0){
      star[isCrossEventHorizon, 1:numVar] <- generateRandom(numStarThatCrossEventHorizon, numVar, lowerBound, upperBound)
      star$fitness <- calcFitness(FUN, optimType, data.matrix(star[,1:numVar]))
    }
    
    # combine blackhole and star 
    candidateSolutions <- rbind(blackhole, star) 
  }
  #??
  return(calcBest(FUN, -1*optimType, as.matrix(candidateSolutions[, 1:numVar])))
}
