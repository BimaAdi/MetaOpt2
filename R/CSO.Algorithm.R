# Cat Swarm Optimization (CSO)

source('./R/metaheuristic.FunctionCollection.R')

CSO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
                mixtureRatio=0.5, tracingConstant=0.1, maximumVelocity=1, smp=3, srd=50, cdc=1, spc=TRUE){
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
  bestPos <- engineCSO(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution,
                       mixtureRatio, tracingConstant, maximumVelocity, smp, srd, cdc, spc)
  return(bestPos)
}

engineCSO <- function(FUN, optimType, numVar, numPopulation, maxIter, lowerBound, upperBound, candidateSolution,
                      mixtureRatio, tracingConstant, maximumVelocity, smp, srd, cdc, spc){
  # generate candidate solutions
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  velocity <- apply(matrix(rep(NA, numPopulation*numVar), ncol = numVar), c(1, 2), function(x){
    runif(1, min = 0, max = maximumVelocity)
  })
  candidateSolutions <- data.frame(candidateSolution, velocity, fitness)
  bestCandidate <- candidateSolutions[order(candidateSolutions$fitness)[1],]
  
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # Give each candidate solutions flag
    candidateSolutions$flag <- flaging(mixtureRatio, numPopulation)
    # Determine/update best and worst candidate solution
    bestCandidateinThisIteration <- candidateSolutions[order(candidateSolutions$fitness)[1],]
    if(bestCandidate$fitness > bestCandidateinThisIteration$fitness) bestCandidate <- bestCandidateinThisIteration
    worstCandidate <- candidateSolutions[order(candidateSolutions$fitness)[numPopulation],]
    indexVariable <- 1:numVar # index column variable in data frame
    indexVelocity <- (numVar+1):(numVar+numVar)# index column velocity in data frame
    
    # if flag == "tracing" ----
    tracing <- candidateSolutions[candidateSolutions$flag == "tracing",]
    tracingVariable <- as.matrix(tracing[,indexVariable])
    tracingVelocity <- as.matrix(tracing[,indexVelocity])
    # Update velocity using
    randomMatrix <- apply(tracingVelocity, c(1, 2), function(x){
      runif(1)
    })
    bestCandidateVariable <- as.numeric(bestCandidate[,1:numVar])
    bestCandidateVariable <- matrix(rep(bestCandidateVariable, nrow(tracing)), ncol = numVar, byrow = TRUE)
    tracingVelocity <- tracingVelocity + randomMatrix * tracingConstant * (bestCandidateVariable - tracingVariable)
    # check velocity
    tracingVelocity[tracingVelocity > maximumVelocity] <- maximumVelocity
    # update position
    tracingVariable <- tracingVariable + tracingVelocity
    # update candidate solution
    candidateSolutions[candidateSolutions$flag == "tracing", indexVariable] <- tracingVariable
    candidateSolutions[candidateSolutions$flag == "tracing", indexVelocity] <- tracingVelocity
    
    
    # if flag == "seeking" ----
    seeking <- candidateSolutions[candidateSolutions$flag == "seeking",]
    seekingVariable <- seeking[,indexVariable]
    if(numVar == 1){
      x <- seekingVariable
      copyId <- 1:length(seekingVariable)
      seekingVariable <- data.frame(x, copyId)
    }else{
      seekingVariable$copyId <- 1:nrow(seekingVariable)
    }
    
    # make smp copies
    copies <- data.frame()
    if(spc == TRUE){
      for(i in 1:smp){
        copies <- rbind(copies, seekingVariable)
      }
    }else{
      for(i in 1:(smp-1)){
        copies <- rbind(copies, seekingVariable)
      }
    }
    
    # modified copies
    if(cdc != 0){
      modified <- apply(as.matrix(copies[,indexVariable]), c(1), function(x){
        pickedVariables <- sample(1:numVar, cdc)
        posOrNeg <- sapply(1:cdc, function(y){
          sample(c(1, -1), 1)
        })
        x[pickedVariables] <- x[pickedVariables]*posOrNeg*srd/100
        return(x)
      })
      if(numVar == 1) copies[,indexVariable] <- modified else copies[,indexVariable] <- t(modified)
    }
    
    # calculate probabilty of all candidate solution (flag == "seeking")
    copies <- rbind(seekingVariable, copies)
    copies$probability <- probability(as.matrix(copies[,indexVariable]), bestCandidate, worstCandidate, FUN, optimType)
    
    # chose one candidate solution for each copy based on probability
    for(i in 1:nrow(seekingVariable)){
      numCopies <- nrow(copies[copies$copyId == i,])
      probCopies <- copies[copies$copyId == i, "probability"]
      choosenCopy <- sample(1:numCopies, 1, prob = probCopies)
      choosenCopy <- copies[copies$copyId == i, ][choosenCopy, indexVariable]
      seekingVariable[seekingVariable$copyId == i, indexVariable] <- choosenCopy
    }
    
    # update candidate solution 
    candidateSolutions[candidateSolutions$flag == "seeking", indexVariable] <- seekingVariable[,indexVariable]
    
    # update candidate solution fitness
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,indexVariable]))
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(as.matrix(bestCandidate[,indexVariable]))
}

flaging <- function(mixtureRatio, numPopulation){
  numSeeking <- mixtureRatio * numPopulation
  numTracing <- (1 - mixtureRatio) * numPopulation
  seeking <- rep("seeking", ceiling(numSeeking))
  tracing <- rep("tracing", ceiling(numTracing))
  result <- sample(c(seeking, tracing), replace = FALSE)
  if(length(result) != numPopulation){
    result <- result[1:numPopulation]
  }
  return(result)
}

probability <- function(input, best, worst, FUN, optimType){
  inputFitness <- calcFitness(FUN, optimType, input)
  bestFitness <- best$fitness
  worstFitness <- worst$fitness
  result <- 1 - abs(inputFitness - bestFitness)/abs(bestFitness - worstFitness)
  return(result)
}

