# Gravitational Based Search Algorithm (GBS)

GBS <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
                gravitationalConst=100, kbest=0.5){
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
  
  # generate initial population
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineGBS(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, 
                       gravitationalConst, kbest)
  return(bestPos)
}

engineGBS <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      gravitationalConst=100, kbest=0.5){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  # give every candidate solution initial velocity
  velocity <- matrix(rep(0, numPopulation * numVar), ncol = numVar)
  
  gbest <- calcBest(FUN, -1*optimType, candidateSolution)
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # update gravitational constant
    eulerConst <- 0.5772156649
    gravitationalConstT <- gravitationalConst / exp(0.01 * t)
    # get best and worst candidate solution in this current t
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1],]
    worst <- candidateSolution[order(CSFitness)[numPopulation],]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))
    
    # calculate gravitaional mass for each candidate solution
    totalGM <- sum((CSFitness - worstFitness) / (bestFitness - worstFitness))
    if(is.nan(totalGM)){
      print("out of bound")
      break
    }
    GM <- (CSFitness - worstFitness) / (bestFitness - worstFitness)
    GM <- GM / totalGM
    
    # calculate total force
    epsilon <- 8.854e-12
    k <- round(numPopulation * kbest)
    candidateSolution <- as.matrix(candidateSolution[as.numeric(order(CSFitness)), ])
    velocity <- as.matrix(velocity[as.numeric(order(CSFitness)), ])
    totalForce <- c()
    for(i in 1:numPopulation){
      totalI <- rep(0, numVar)
      for(j in 1:k){
        distance <- as.numeric(dist(rbind(candidateSolution[i,],candidateSolution[j,])))
        Force <- runif(1) * gravitationalConstT * GM[j] * (candidateSolution[j,] - candidateSolution[i, ])/distance
        if(all(is.nan(Force))){
          Force <- as.numeric(rep(0, numVar))
        }
        totalI <- totalI + Force
      }
      totalForce <- rbind(totalForce, totalI)
    }
    
    randomMatrix <- matrix(runif(numPopulation * numVar), ncol = numVar)
    velocity <- randomMatrix * velocity + totalForce
    candidateSolution <- candidateSolution + velocity
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(gbest)  
}
