# Krill-Heard Algorithm(KH)

source('./R/metaheuristic.FunctionCollection.R')

FUN <- sphere
optimType <- "MIN"
numVar <- 1
numPopulation <- 10
maxIter <- 3
rangeVar <- matrix(c(-10, 10), nrow = 2)
maxMotionInduced <- 0.01
inertiaWeightOfMotionInduced <- 0.5
epsilon <- 5e-05
foragingSpeed <- 0.02
inertiaWeightOfForagingSpeed <- 0.2 #[0, 1]
maxDifussionSpeed <- 0.001 # [0.002, 0.01]
constantSpace <- 0.1 #[0, 2]
mu <- 0.5 #[0, 1]

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

# generate candidate solutions
indexVariable <- 1:numVar
fitness <- calcFitness(FUN, optimType, candidateSolution)
motionInduced <- matrix(rep(0, numPopulation * numVar), ncol = numVar)
foragingMotion <- matrix(rep(0, numPopulation * numVar), ncol = numVar)
physicalDifusion <- matrix(rep(0, numPopulation * numVar), ncol = numVar) 
candidateSolutions <- data.frame(candidateSolution, fitness)
candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,1:numVar]))
best <- candidateSolutions[order(candidateSolutions$fitness)[1],]
worst <- candidateSolutions[order(candidateSolutions$fitness)[numPopulation],]

for(t in 1:maxIter){
  print(t)
  # Motion Induced Calculation ----
  # alpha local
  sensingDistance <- isIncludeSensingDistances(candidateSolutions)
  alphaLocal <- c()
  for(index in 1:numPopulation){
    # motion calculation
    jVariable <- as.matrix(candidateSolutions[sensingDistance[index,], indexVariable])
    iVariable <- as.numeric(candidateSolutions[index, indexVariable])
    motion <- KHMotion(iVariable, jVariable, epsilon)
    
    # fitness calculation
    jFitness <- as.matrix(candidateSolutions[sensingDistance[index, ], "fitness"])
    iFitness <- candidateSolutions[index, "fitness"]
    fitness <- KHFitnessMotion(iFitness, jFitness, worst$fitness, best$fitness)
    
    # alpha local calculation
    alphaLocal <- rbind(alphaLocal, colSums(fitness * motion))
  }
  
  random <- randomMatrix(numPopulation, numVar, 0, 1)
  # as.numeric ??
  KHMotionCSBest <- KHMotion(as.matrix(candidateSolutions[,indexVariable]), as.numeric(best[,indexVariable]), epsilon)
  KHFitnessMotionCSBest <- KHFitnessMotion(candidateSolutions[,"fitness"], best$fitness, worst$fitness, best$fitness)
  motionInduced <- maxMotionInduced * (alphaLocal + 2 * (random + t/maxIter) * KHFitnessMotionCSBest * KHMotionCSBest) + 
    inertiaWeightOfMotionInduced * motionInduced
  
  # foraging motion calculation ----
  if(numVar == 1){
    food <- sum(1 / candidateSolutions$fitness * candidateSolutions[,indexVariable])/sum(1 / candidateSolutions$fitness)
  }else{
    food <- colSums(1 / candidateSolutions$fitness * candidateSolutions[,indexVariable])/sum(1 / candidateSolutions$fitness)
  }
  foodFitness <- calcFitness(FUN, optimType, t(food))
  betaFood <- KHFitnessMotion(candidateSolutions[,"fitness"], foodFitness, worst$fitness, best$fitness) * 
    KHMotion(as.matrix(candidateSolutions[,indexVariable]), food, epsilon) 
  betaBest <- KHFitnessMotion(candidateSolutions[,"fitness"], best$fitness, worst$fitness, best$fitness) *
    KHMotion(as.matrix(candidateSolutions[,indexVariable]), best[,indexVariable], epsilon)
  foragingMotion <- foragingSpeed * (2*(1 - t/maxIter) * betaFood * betaBest) + inertiaWeightOfForagingSpeed * foragingMotion
  
  # Physical Difusion Calculation ----
  randomDirectionalVector <- randomMatrix(numPopulation, numVar, -1, 1)
  physicalDifusion[,] <- maxDifussionSpeed * (1 - t/maxIter) * randomDirectionalVector
  
  # calculate motion calculation & update candidate solutions ----
  motionCalculation <- motionInduced + foragingMotion + physicalDifusion
  deltaT <- constantSpace * sum(upperBound - lowerBound)
  candidateSolutions[,indexVariable] <- candidateSolutions[,indexVariable] + deltaT * motionCalculation
  candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,1:numVar]))
  best <- candidateSolutions[order(candidateSolutions$fitness)[1],]
  worst <- candidateSolutions[order(candidateSolutions$fitness)[numPopulation],]
  
  # implemented genetic operator crossover
  idBest <- order(candidateSolutions$fitness)[1]
  print(KHFitnessMotion(candidateSolutions[-1*idBest,"fitness"], best$fitness, worst$fitness, best$fitness))
  Cr <- 0.2 * KHFitnessMotion(candidateSolutions[-1*idBest,"fitness"], best$fitness, worst$fitness, best$fitness)
  print(Cr)
  Cr <- matrix(as.numeric(rep(Cr, numVar)), ncol = numVar)
  random <- randomMatrix(numPopulation - 1, numVar, 0, 1)
  probCr <- random < Cr
  if(!all(probCr == FALSE)){
    print("in probCr")
    CrCandidateSolutions <- as.matrix(candidateSolutions[-1*idBest, indexVariable])
    CrCandidateSolutions[probCr] <- sapply(col(probCr)[probCr], function(x){
      choosenCS <- sample(1:(numPopulation - 1), 1)
      return(CrCandidateSolutions[choosenCS, x])
    })
    candidateSolutions[-1*idBest, indexVariable] <- CrCandidateSolutions
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,1:numVar]))
    best <- candidateSolutions[order(candidateSolutions$fitness)[1],]
    worst <- candidateSolutions[order(candidateSolutions$fitness)[numPopulation],]
  }
  
  # implemented genetic operator mutation
  idBest <- order(candidateSolutions$fitness)[1]
  Mu <- 0.05/KHFitnessMotion(candidateSolutions[-1*idBest,"fitness"], best$fitness, worst$fitness, best$fitness)
  Mu <- matrix(as.numeric(rep(Mu, numVar)), ncol = numVar)
  random <- randomMatrix(numPopulation - 1, numVar, 0, 1)
  probMu <- random < Mu
  if(!all(probMu == FALSE)){
    MuCandidateSolutions <- as.matrix(candidateSolutions[-1*idBest, indexVariable])
    MuCandidateSolutions[probMu] <- sapply(col(probMu)[probMu], function(x){
      bestVariable <- best[,indexVariable][,x]
      item <- 1:(numPopulation - 1)
      p <- MuCandidateSolutions[sample(item, 1), x]
      q <- MuCandidateSolutions[sample(item, 1), x]
      return(bestVariable + mu*(p - q))
    })
    candidateSolutions[-1*idBest, indexVariable] <- MuCandidateSolutions
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,1:numVar]))
    best <- candidateSolutions[order(candidateSolutions$fitness)[1],]
    worst <- candidateSolutions[order(candidateSolutions$fitness)[numPopulation],]
  }
}


isIncludeSensingDistances <- function(input){
  test <- as.matrix(input)
  distances <- as.matrix(dist(test, diag = TRUE, upper = TRUE))
  n <- ncol(distances)
  sensingDistances <- apply(distances, c(1), function(x){
    1/5*n*sum(x)
  })
  return(distances < sensingDistances)
}

# KHMotion
KHMotion <- function(i, j, epsilon){
  # either or both i and j must be matrix with same size
  if(!is.matrix(i) & !is.matrix(j)){
    stop("either or both i and j must be matrix with same size")
  }else if(!is.matrix(i)){
    i <- matrix(rep(as.numeric(i), nrow(j)), ncol = ncol(j), byrow = TRUE)
  }else if(!is.matrix(j)){
    j <- matrix(rep(as.numeric(j), nrow(i)), ncol = ncol(i), byrow = TRUE)
  }
  nv <- ncol(i)
  result <- as.matrix(apply(cbind(i, j), c(1), function(x){
    indexI <- 1:nv
    indexJ <- (nv+1):(nv*2)
    return(x[indexJ] - x[indexI]/dist(rbind(x[indexJ], x[indexI])) + epsilon)
  }))
  if(nv > 1){
    result <- t(result)
  }
  result[result == Inf | result == -Inf] <- 0
  return(result)
}

# KHFitnessMotion
KHFitnessMotion <- function(i, j, worst, best){
  # fitness
  if(length(i) == 1 & length(j) > 1){
    i <- rep(as.numeric(i), length(j))
  }else if(length(i) > 1 & length(j) == 1){
    j <- rep(as.numeric(j), length(i))
  }
  result <- as.numeric(i - j/worst - best)
  result[is.nan(result)] <- -Inf
  return(result)
}

randomMatrix <- function(numberRow, numberColumn, minRange, maxRange){
  apply(matrix(1:(numberColumn*numberRow), ncol = numberColumn, nrow = numberRow), c(1, 2), function(x){
    runif(1, min = minRange, max = maxRange)
  })
}
