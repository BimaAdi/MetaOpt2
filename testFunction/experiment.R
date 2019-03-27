FUN <- f1
optimType <- "MIN"
numVar <- 2
numPopulation <- 10
maxIter <- 10
rangeVar <- matrix(c(-10, 10), nrow = 2)
selectionSize=10
multipicationFactor=0.5
hypermutationRate=0.5

lowerBound <- c(-10, -20)
upperBound <- c(10, 20)
x <- matrix(1:10, ncol = 2, byrow = TRUE)
probMatrix <- matrix(NA, ncol = 2, nrow= 5, byrow = TRUE)
probMatrix <- apply(matrix(NA, ncol = 2, nrow= 5, byrow = TRUE), c(1, 2), function(x){
  runif(1)
})
x[probMatrix < 0.5]   
col(x)[probMatrix < 0.5]

lowerBound[col(x)[probMatrix < 0.5]]
upperBound[col(x)[probMatrix < 0.5]]
mutationMatrix <- apply(rbind(lowerBound[col(x)[probMatrix < 0.5]], upperBound[col(x)[probMatrix < 0.5]]), c(2), function(x){
  runif(1, min = x[1], max = x[2])
})
x[probMatrix < 0.5] <- mutationMatrix

# 1. pseudocode diganti atau tidak
# 2. kecepatan algoritma DE
# 3. apakah format program sudah sesuai

z <- function(x){
  return(sum(x^2))
}

metaOpt(z, optimType = "MIN", algorithm = c("ABC", "BHO"), numVar = 5, rangeVar = matrix(c(-100, 100), nrow=2), seed = 1)
BHO(FUN = z, optimType = "MIN", numVar=2, rangeVar = matrix(c(-100, 100), nrow=2))

CS <- generateRandom(3, 2, -10, 10)
CS <- matrix(c(1, -1, -1, 0, 0, 1), ncol = 2, byrow = TRUE)
lb <- rep(-10, 2)
ub <- rep(10, 2)
engineBHO(FUN = z, optimType = 1, maxIter = 500, lowerBound = lb, upperBound = ub, candidateSolution = CS)

set.seed(2)
SFL(FUN=f2, optimType = "MIN", numVar = 10, numPopulation = 50, maxIter = 500, rangeVar = matrix(c(-10, 10), nrow = 2))
set.seed(2)
BHO(FUN = f2, optimType = "MIN", numVar = 10, numPopulation = 50, maxIter = 500, rangeVar = matrix(c(-10, 10), nrow = 2))

calcFitness(f12, 1, t(test1))
calcFitness(f12, 1, t(test2))   

set.seed(20)
FUN <- f1
rangeVar <- matrix(c(-100, 100), nrow = 2)
kh1 <- KH(FUN = FUN, optimType = "MIN", numVar = 10, numPopulation = 50, maxIter = 100, rangeVar = rangeVar)
kh2 <- KH(FUN = FUN, optimType = "MIN", numVar = 10, numPopulation = 50, maxIter = 100, rangeVar = rangeVar, maxMotionInduced = 0.01, foragingSpeed = 0.1, maxDifussionSpeed = 0.1)
kh3 <- KH(FUN = FUN, optimType = "MIN", numVar = 10, numPopulation = 50, maxIter = 100, rangeVar = rangeVar, maxMotionInduced = 0.01, foragingSpeed = 0.01, maxDifussionSpeed = 0.01)
calcFitness(FUN, 1, t(kh1))
calcFitness(FUN, 1, t(kh2))
calcFitness(FUN, 1, t(kh3))


tic()
SFL(FUN = FUN, optimType = "MIN", numVar = 10, numPopulation = 10, maxIter = 1000, rangeVar = rangeVar)
toc()
