# Krill-Heard Algorithm(KH)

KH <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
               maxMotionInduced=0.01, inertiaWeightOfMotionInduced=0.01, epsilon=1e-05, foragingSpeed=0.02,
               inertiaWeightOfForagingSpeed=0.01, maxDifussionSpeed=0.01, constantSpace=1, mu=0.1){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }
  
  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }
  
  if(inertiaWeightOfMotionInduced < 0 | inertiaWeightOfMotionInduced > 1){
    stop("inertiaWeightOfMotionInduced must between 0 and 1")
  }
  
  if(inertiaWeightOfForagingSpeed < 0 | inertiaWeightOfForagingSpeed > 1){
    stop("inertiaWeightOfForagingSpeed must between 0 and 1")
  }
  
  if(maxMotionInduced < 0 | maxMotionInduced > 1){
    stop("maxMotionInduced must between 0 and 1")
  }
  
  if(foragingSpeed < 0 | foragingSpeed > 1){
    stop("foragingSpeed must between 0 and 1")
  }
  
  if(maxDifussionSpeed < 0 | maxDifussionSpeed > 1){
    stop("maxDifussionSpeed must between 0 and 1")
  }
  
  if(constantSpace < 0 | constantSpace > 2){
    stop("constantSpace must between 0 and 2")
  }
  
  if(mu < 0 | mu > 1){
    stop("mu must between 0 and 1")
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
  bestPos <- engineKH(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      maxMotionInduced, inertiaWeightOfMotionInduced, epsilon, foragingSpeed,
                      inertiaWeightOfForagingSpeed, maxDifussionSpeed, constantSpace, mu)
  return(bestPos)
}

engineKH <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                     maxMotionInduced, inertiaWeightOfMotionInduced, epsilon, foragingSpeed,
                     inertiaWeightOfForagingSpeed, maxDifussionSpeed, constantSpace, mu){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  
  N <- matrix(rep(0, numPopulation * numVar), ncol = numVar) 
  f <- matrix(rep(0, numPopulation * numVar), ncol = numVar)   
  
  gbest <- calcBest(FUN, -1*optimType, candidateSolution)
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1], ]
    worst <- candidateSolution[order(CSFitness)[numPopulation], ]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
    
    
    # motion iduced
    sensingDistance <- 1/5*ncol(candidateSolution)*colSums(as.matrix(dist(candidateSolution)))
    isIncludeSD <- as.matrix(dist(candidateSolution, diag = T, upper = T)) < sensingDistance
    
    alpha <- c()
    for(index in 1:numPopulation){
      X <- apply(as.matrix(candidateSolution[isIncludeSD[index,],]), c(1), function(x, y){
        xijKH(y, x, epsilon)
      }, y=candidateSolution[index,])
      
      K <- sapply(CSFitness[isIncludeSD[index,]], function(x, y){
        kijKH(y,x, bestFitness, worstFitness)
      }, y=CSFitness[index])
      
      if(numVar == 1){
        alphaLocal <- sum(X * K)
      }else{
        alphaLocal <- colSums(t(X) * K)
      }
      
      Cbest <- 2*(runif(1)+t/maxIter)
      X <- xijKH(candidateSolution[index,], best, epsilon)
      K <- kijKH(CSFitness[index], bestFitness, bestFitness, worstFitness)
      alphaTarget <- Cbest*K*X
      alpha <- rbind(alpha, alphaLocal + alphaTarget)  
    }
    
    N <- maxMotionInduced * alpha + inertiaWeightOfMotionInduced * N
    
    # foraging motion
    if(numVar == 1){
      Xfood  <- sum(candidateSolution * 1 / CSFitness) / sum(1/CSFitness)
      if(is.nan(Xfood)){
        Xfood <- 0
      }
    }else{
      Xfood  <- colSums(candidateSolution * 1 / CSFitness) / sum(1/CSFitness)
    }
    xFoodFitness <- calcFitness(FUN, optimType, matrix(Xfood, ncol = numVar))
    Cfood <- 2*(1-t/maxIter)
    Kifood <- sapply(CSFitness, function(x){
      kijKH(x, xFoodFitness, bestFitness, worstFitness)  
    })
    Xifood <- apply(candidateSolution, c(1), function(x){
      xijKH(x, Xfood, epsilon)
    })
    
    Kibest <- sapply(CSFitness, function(x){
      kijKH(x, bestFitness, bestFitness, worstFitness)  
    })
    Xibest <- apply(candidateSolution, c(1), function(x){
      xijKH(x, best, epsilon)
    })
    
    if(numVar == 1){
      betaFood <- Cfood*Kifood*Xifood
      betaBest <- Xibest * Kibest
    }else{
      betaFood <- t(Cfood*Kifood*Xifood)
      betaBest <- t(Xibest) * Kibest
    }
    beta <- betaFood + betaBest
    
    f <- foragingSpeed*beta + inertiaWeightOfForagingSpeed*f
    
    # physical difussion
    D <- maxDifussionSpeed * (1 - t/maxIter)*runif(1, min = -1, max = 1)
    
    # Motion calculation
    TotalMotion <- N + f + D
    deltaT <- constantSpace*sum(upperBound - lowerBound)
    candidateSolution <- candidateSolution + deltaT * TotalMotion
    
    # implement genetic operator ----
    
    # CrossOver
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1], ]
    worst <- candidateSolution[order(CSFitness)[numPopulation], ]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
    Kibest <- sapply(CSFitness, function(x){
      kijKH(x, bestFitness, bestFitness, worstFitness)  
    })
    randomMatrix <- matrix(runif(numVar * numPopulation), ncol = numVar)
    Cr <- 0.2 * Kibest
    prob <- randomMatrix < Cr
    if(!all(prob == FALSE)){
      Xrm <- sapply(col(candidateSolution)[prob], function(x){
        choosen <- sample.int(numPopulation, 1)
        return(candidateSolution[choosen, x])
      })
      candidateSolution[prob] <- Xrm
    }
    
    # Mutation
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1], ]
    worst <- candidateSolution[order(CSFitness)[numPopulation], ]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
    Kibest <- sapply(CSFitness, function(x){
      kijKH(x, bestFitness, bestFitness, worstFitness)  
    })
    randomMatrix <- matrix(runif(numVar * numPopulation), ncol = numVar)
    Mu <- 0.05 * Kibest
    prob <- randomMatrix < Mu
    if(!all(prob == FALSE)){
      Xgbest <- sapply(col(candidateSolution)[prob], function(x){
        P <- sample.int(numPopulation, 1)
        Q <- sample.int(numPopulation, 1)
        return(gbest[x] + mu * (candidateSolution[P, x] - candidateSolution[Q, x]))
      })
      candidateSolution[prob] <- Xgbest
    }
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
  return(gbest)
}


xijKH <- function(i, j, epsilon){
  (j - i)/(dist(rbind(j, i)) + epsilon)
}

kijKH <- function(i, j, best, worst){
  if(worst == best){
    0  
  }else{
    (i - j)/(worst - best)
  }
}
