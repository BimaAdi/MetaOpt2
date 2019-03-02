# Cuckoo Search (CS)

CS <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, abandonedFraction=0.5){
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
  
  #generate candidate solutions
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineCS(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, abandonedFraction)
  
  return(bestPos)
}

engineCS <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, abandonedFraction){
  numPopulation <- nrow(candidateSolution)
  numVar <- ncol(candidateSolution)
  best <- c()
  
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
    # create a cuckoo egg
    choosen <- sample.int(numPopulation, 1)
    choosen <- candidateSolution[choosen,]
    cuckooegg <- choosen + levyFlight(choosen, best)
    
    # compare with candidate solution
    otherChoosen <- sample.int(numPopulation, 1)
    otherFitness <- calcFitness(FUN, optimType, matrix(candidateSolution[otherChoosen, ], ncol = numVar))
    cuckooFitness <- calcFitness(FUN, optimType, matrix(cuckooegg, ncol = numVar))
    if(cuckooFitness < otherFitness){
      candidateSolution[otherChoosen, ] <- cuckooegg
    }
    
    # remove fraction
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    candidateSolution <- as.matrix(candidateSolution[order(CSFitness), ])
    numCSToRemove <- round(numPopulation * abandonedFraction)
    removeIndex <- as.numeric((numPopulation - numCSToRemove + 1):numPopulation)
    candidateSolution[removeIndex,] <- generateRandom(numCSToRemove, numVar, lowerBound, upperBound)
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
  return(best)  
}

levyFlight <- function(CS, best){
  n <- length(CS)
  beta <- 3/2
  sigma <- (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta)
  u <- rnorm(n) * sigma
  v <- rnorm(n)
  step <- u/abs(v)^(1/beta)
  stepsize <- 0.01 * step * (CS - best)
  stepsize <- stepsize * rnorm(n) 
  return(stepsize)
}
