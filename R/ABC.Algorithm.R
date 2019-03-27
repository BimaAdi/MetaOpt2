# Artificial Bee Colony Algorithm (ABC)

ABC <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, cycleLimit=as.integer(10)){
  
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }
  
  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }
  
  if(cycleLimit < 0){
    stop("cycleLimit must greater than 0")
  }else if(!is.integer(cycleLimit)){
    stop("cycleLimit must be integer (as.integer())")
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
  bestPos <- engineABC(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, cycleLimit)
  
  return(bestPos)
}

engineABC <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, cycleLimit){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  # generate candidate solutions
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  limit <- rep(cycleLimit, numPopulation)
  candidateSolutions <- data.frame(candidateSolution, fitness, limit)
  
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # employeed bee phase
    for(i in 1:numPopulation){
      pickedColumn <- sample(1:numVar, 1)
      temp <- 1:numPopulation
      temp <- temp[!temp %in% c(i)]
      pickedRow <- sample(temp, 1)
      new <- as.matrix(candidateSolutions[i, 1:numVar])
      new[1, pickedColumn] <- candidateSolutions[i, pickedColumn] + 
        runif(1, min=-1, max=1)*(candidateSolutions[i, pickedColumn] - candidateSolutions[pickedRow, pickedColumn])
      if(calcFitness(FUN, optimType, new) < calcFitness(FUN, optimType, as.matrix(candidateSolutions[i, 1:numVar]))){
        candidateSolutions[i, 1:numVar] <- new
        candidateSolutions[i, "limit"] <- cycleLimit
      }else{
        candidateSolutions[i, "limit"] <- candidateSolutions[i, "limit"] - 1
      }
    }
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,1:numVar]))
    # prevent best candidate solution abadnoned
    if(candidateSolutions[order(candidateSolutions$fitness)[1],"limit"] <= 0){
      candidateSolutions[order(candidateSolutions$fitness)[1],"limit"] <- cycleLimit
    }
    
    # onlooker bee phase
    for(i in 1:numPopulation){
      if(runif(1) < (candidateSolutions$fitness[i]/sum(candidateSolutions$fitness))){
        pickedColumn <- sample(1:numVar, 1)
        temp <- 1:numPopulation
        temp <- temp[!temp %in% c(i)]
        pickedRow <- sample(temp, 1)
        new <- as.matrix(candidateSolutions[i, 1:numVar])
        new[1, pickedColumn] <- candidateSolutions[i, pickedColumn] + 
          runif(1, min=-1, max=1)*(candidateSolutions[i, pickedColumn] - candidateSolutions[pickedRow, pickedColumn])
        if(calcFitness(FUN, optimType, new) < calcFitness(FUN, optimType, as.matrix(candidateSolutions[i, 1:numVar]))){
          candidateSolutions[i, 1:numVar] <- new
          candidateSolutions[i, "limit"] <- cycleLimit
        }else{
          candidateSolutions[i, "limit"] <- candidateSolutions[i, "limit"] - 1
        }
      }
    }
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,1:numVar]))
    # prevent best candidate solution abadnoned
    if(candidateSolutions[order(candidateSolutions$fitness)[1],"limit"] <= 0){
      candidateSolutions[order(candidateSolutions$fitness)[1],"limit"] <- cycleLimit
    }
    
    # scout bee phase
    # check abandoned candidate solution
    numAbandoned <- nrow(candidateSolutions[candidateSolutions$limit <= 0,])
    if(numAbandoned != 0){
      candidateSolution <- generateRandom(numAbandoned, numVar, lowerBound, upperBound)
      fitness <- calcFitness(FUN, optimType, candidateSolution)
      limit <- rep(cycleLimit, numAbandoned)
      candidateSolutions[candidateSolutions$limit <= 0,] <- data.frame(candidateSolution, fitness, limit)
    }
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, -1*optimType, as.matrix(candidateSolutions[, 1:numVar])))
}
