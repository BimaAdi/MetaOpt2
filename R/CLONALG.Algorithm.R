# Clonal Selection Algorithm (CLONALG)

CLONALG <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, 
                    selectionSize=10, multipicationFactor=0.5, hypermutationRate=0.5){
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
  candidateSolutions <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineCLONALG(FUN, optimType, maxIter, rangeVar, lowerBound, upperBound,
                           selectionSize, multipicationFactor, hypermutationRate, candidateSolutions)
  return(bestPos)
}

engineCLONALG <- function(FUN, optimType, maxIter, rangeVar, lowerBound, upperBound,
                          selectionSize, multipicationFactor, hypermutationRate,
                          candidateSolution){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  # evaluate candidate solution
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  candidateSolutions <- data.frame(candidateSolution, fitness)
  
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # Select top "selectionSize" with best fitness from candidateSolutions as topSelections
    candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness), ]
    topSelections <- data.matrix(candidateSolutions[1:selectionSize, 1:numVar])
    
    # clone topSelections as clone
    # create empty matrix
    clone <- matrix()[-1:-1]
    for(i in 1:selectionSize){
      numClone <- round(multipicationFactor*numPopulation/i, digits = 0)
      clone <- rbind(clone, matrix(data=rep(topSelections[i,], numClone), ncol = numVar, byrow = TRUE))
    }
    
    # hypermutate clone
    probMatrix <- apply(matrix(NA, ncol = ncol(clone), nrow= nrow(clone), byrow = TRUE), c(1, 2), function(x){
      runif(1)
    })
    rangeVarMatrix <- rbind(lowerBound[col(clone)[probMatrix <= multipicationFactor]], 
                            upperBound[col(clone)[probMatrix <= multipicationFactor]])
    clone[probMatrix <= multipicationFactor] <- apply(rangeVarMatrix, c(2), function(x){
      runif(1, min = x[1], max = x[2])
    })
    
    # maturation step
    candidateSolution <- clone
    fitness <- calcFitness(FUN, optimType, candidateSolution)
    candidateSolutions <- rbind(candidateSolutions, data.frame(candidateSolution, fitness))
    candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness), ]
    candidateSolutions <- candidateSolutions[1:numPopulation, ]
    
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, optimType, as.matrix(candidateSolutions[, 1:numVar])))
}
