# Cuckoo Search (CS)
require("rmutil")
source('./R/metaheuristic.FunctionCollection.R')

FUN <- sphere
optimType <- "MIN"
numPopulation <- 25
numVar <- 1
maxIter <- 500
rangeVar <- matrix(c(-10, 10), nrow = 2)
stepSize <- 1
abandonedFraction <- 0.5

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

fitness <- calcFitness(FUN, optimType, candidateSolution)
candidateSolutions <- data.frame(candidateSolution, fitness)
indexVariable <- 1:numVar
numRemove <- round(numPopulation*abandonedFraction)

for(t in 1:maxIter){
  #print(t)
  # move a candidate solution randomly as cuckoo egg
  choosenCS <- sample(1:numPopulation, 1)
  choosenCS <- as.matrix(candidateSolutions[choosenCS, indexVariable])
  levyRandom <- rlevy(n=numVar)
  levyRandom[levyRandom > 3] <- 3
  levyRandom <- matrix(levyRandom, ncol = numVar)
  cuckoo <- choosenCS + stepSize * levyRandom
  
  # pick a candidate solution randomly and compare it with cuckoo egg
  pickedCSIndex <- sample(1:numPopulation, 1)
  if(candidateSolutions[pickedCSIndex, "fitness"] > calcFitness(FUN, optimType, as.matrix(cuckoo))){
    candidateSolutions[pickedCSIndex, indexVariable] <- as.matrix(cuckoo)
    candidateSolutions[pickedCSIndex, "fitness"] <- calcFitness(FUN, optimType, as.matrix(cuckoo))
  }
  
  # remove fraction (abandonedFraction) from candidate solutions with worst fitness
  candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness),]
  removeCSIndex <- (numPopulation - numRemove + 1):numPopulation
  bestVar <- as.numeric(candidateSolutions[order(candidateSolutions$fitness)[1],indexVariable])
  candidateSolutions[removeCSIndex, indexVariable] <- generateRandomCS(numRemove, dimension, bestVar, stepSize, lowerBound, upperBound)
  candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,indexVariable]))
  candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness),]
}

generateRandomCS <- function(numRemove, numVar, bestVar, stepSize, lowerBound, upperBound){
  numLocalSearch <- ceiling(numVar/2)
  numPureRandom <- numRemove - numLocalSearch
  # localSearch
  levyRandom <- rlevy(n=numVar*numLocalSearch)
  levyRandom[levyRandom > 3] <- 3
  levyRandom <- matrix(levyRandom, ncol = numVar)
  localSearch <- t(bestVar + stepSize * t(levyRandom))
  # pureRandom
  pureRandom <- generateRandom(numPureRandom, numVar, lowerBound, upperBound)
  result <- as.matrix(rbind(localSearch, pureRandom))
  return(result)
}
# -------------------------------------------------------------------
progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
for(t in 1:maxIter){
  for(index in 1:numPopulation){
    # move a candidate solution randomly as cuckoo egg
    choosenCS <- as.matrix(candidateSolutions[index, indexVariable])
    levyRandom <- rlevy(n=numVar)
    levyRandom[levyRandom > 3] <- 3
    levyRandom <- matrix(levyRandom, ncol = numVar)
    cuckoo <- choosenCS + stepSize * levyRandom
    
    # pick a candidate solution randomly and compare it with cuckoo egg
    pickedCSIndex <- sample(1:numPopulation, 1)
    if(candidateSolutions[pickedCSIndex, "fitness"] > calcFitness(FUN, optimType, as.matrix(cuckoo))){
      candidateSolutions[pickedCSIndex, indexVariable] <- as.matrix(cuckoo)
      candidateSolutions[pickedCSIndex, "fitness"] <- calcFitness(FUN, optimType, as.matrix(cuckoo))
    }
    
    # remove fraction (abandonedFraction) from candidate solutions with worst fitness
    candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness),]
    removeCSIndex <- (numPopulation - numRemove + 1):numPopulation
    candidateSolutions[removeCSIndex, indexVariable] <- generateRandom(numRemove, dimension, lowerBound, upperBound)
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,indexVariable]))
  }
  setTxtProgressBar(progressbar, t)
}
close(progressbar)

