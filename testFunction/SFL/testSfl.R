# Test for SFL------------------------------------------------------------------------
library("tictoc", lib.loc="~/R/win-library/3.5")
source("./testFunction/13TestFunction.R")
source("./R/SFL.Algorithm.R")

inputParameter <- read.csv("./testFunction/SFL/inputSflLampiran.csv", sep = ";")
inputParameter$fungsi <- as.character(inputParameter$fungsi)

test.sfl <- function(inputParameter){
  testResults <- data.frame()
  for(index in 1:nrow(inputParameter)){
    print(cbind(index, nrow(inputParameter)))
    sprintf("%i of %i", index, nrow(inputParameter))
    # prepare parameter
    optimType <- "MIN"
    numVar <- inputParameter$numVar[index]
    rangeVar <- switch(as.character(inputParameter$fungsi[index]), 
                       "f1"= matrix(c(-100, 100), nrow = 2), "f2"=matrix(c(-10, 10), nrow = 2), 
                       "f3"=matrix(c(-100, 100), nrow = 2), "f4"=matrix(c(-100, 100), nrow = 2), 
                       "f5"=matrix(c(-30, 30), nrow = 2), "f6"=matrix(c(-100, 100), nrow = 2), 
                       "f7"=matrix(c(-1.28, 1.28), nrow = 2), "f8"=matrix(c(-500, 500), nrow = 2), 
                       "f9"=matrix(c(-5.12, 5.12), nrow = 2), "f10"=matrix(c(-32, 32), nrow = 2),
                       "f11"=matrix(c(-600, 600), nrow = 2), "f12"=matrix(c(-50, 50), nrow = 2), 
                       "f13"=matrix(c(-50, 50), nrow = 2))
    fungsi <- switch(as.character(inputParameter$fungsi[index]), 
                     "f1"=f1, "f2"=f2, "f3"=f3, "f4"=f4, "f5"=f5, 
                     "f6"=f6, "f7"=f7, "f8"=f8, "f9"=f9, "f10"=f10,
                     "f11"=f11, "f12"=f12, "f13"=f13)
    iterasi <- inputParameter$iterasi[index]
    populasi <- inputParameter$populasi[index]
    numMemeplex <- inputParameter$numMemeplex[index]
    frogLeapingIteration <- inputParameter$frogLeapingIteration[index]
    
    # testing
    tic("finished in")
    x <- SFL(fungsi, optimType = optimType, numVar = numVar, rangeVar = rangeVar,
                 maxIter = iterasi, numPopulation = populasi,
                 numMemeplex = numMemeplex, frogLeapingIteration = frogLeapingIteration)
    time <- toc()
    time <- time$toc - time$tic
    optimum <- calcFitness(fungsi, 1, t(x))
    
    # save result
    x <- matrix(x, nrow = 1)
    if(ncol(x) < 10){
      x <- cbind(x, matrix(rep(0, 10 - ncol(x)), nrow = 1))
    }
    fungsi <- inputParameter$fungsi[index]
    testResult <- data.frame(fungsi, iterasi, populasi, numVar,
                                    numMemeplex, frogLeapingIteration,
                                    optimum, x, time)
    testResults <- rbind(testResults, testResult)
  }
  rownames(testResults) <- 1:nrow(testResults)
  return(testResults)
}

result <- test.sfl(inputParameter)
write.csv2(result, file = "./testFunction/SFL/resultSflLampiran.csv")
