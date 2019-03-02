# Bat Algorithm (BA)

BA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, 
               maxFrequency=1, minFrequency=-1, gama=1, alpha=0.5){
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
  bestPos <- engineBA(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      maxFrequency, minFrequency, gama, alpha)
  
  return(bestPos)
}

engineBA <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                     maxFrequency, minFrequency, gama, alpha){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  
  pf <- matrix(runif(numPopulation * numVar), ncol = numVar)
  velocity <- matrix(rep(0, numPopulation * numVar), ncol = numVar)
  A <- runif(numPopulation)
  pr <- runif(numPopulation)
  
  best <- c()
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    
    # update pf, velocity and candidate solution
    randomMatrix <- matrix(runif(numPopulation * numVar), ncol = numVar)
    pf <- minFrequency + randomMatrix *(maxFrequency - minFrequency)
    velocity <- velocity + t(t(candidateSolution) - best) * pf
    candidateSolution <- candidateSolution + velocity
    candidateSolution <- checkBound(candidateSolution, lowerBound, upperBound)
    
    # pr comparison
    prob <- runif(numPopulation) > pr
    if(!all(prob == FALSE)){
      candidateSolution[prob, ] <- outer(A[prob], best, FUN = "*")
    }
    
    # loudness comparison phase
    fitness <- calcFitness(FUN, optimType, candidateSolution)
    prob <- runif(numPopulation) < A & fitness < bestFitness
    if(!all(prob == FALSE)){
      candidateSolution[prob, ] <- generateRandom(length(prob[prob == TRUE]), numVar, lowerBound, upperBound)
      pr[prob] <- pr[prob]*(1 - exp(-1*gama))
      A[prob] <- alpha*A[prob]
    }
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, -1*optimType, rbind(candidateSolution, best)))
}
