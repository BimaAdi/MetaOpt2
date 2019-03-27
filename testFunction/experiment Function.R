source("./R/metaheuristic.FunctionCollection.R")
source("./R/metaheuristic.mainFunction.R")
source("./testFunction/13TestFunction.R")

experiment <- function(metaAlgorithm, inputFile, outputFile){
  optimType <- "MIN"
  numVar <- 10
  input <- read.csv2(inputFile)
  output <- data.frame()
  
  for(index in 1:nrow(input)){
    print(cbind(index, nrow(input)))
    # preparing parameter
    fungsi <- switch(as.character(input$fungsi[index]), 
                     "f1"=f1, "f2"=f2, "f3"=f3, "f4"=f4, "f5"=f5, 
                     "f6"=f6, "f7"=f7, "f8"=f8, "f9"=f9, "f10"=f10,
                     "f11"=f11, "f12"=f12, "f13"=f13)
    rangeVar <- switch(as.character(input$fungsi[index]), 
                       "f1"= matrix(c(-100, 100), nrow = 2), "f2"=matrix(c(-10, 10), nrow = 2), 
                       "f3"=matrix(c(-100, 100), nrow = 2), "f4"=matrix(c(-100, 100), nrow = 2), 
                       "f5"=matrix(c(-30, 30), nrow = 2), "f6"=matrix(c(-100, 100), nrow = 2), 
                       "f7"=matrix(c(-1.28, 1.28), nrow = 2), "f8"=matrix(c(-500, 500), nrow = 2), 
                       "f9"=matrix(c(-5.12, 5.12), nrow = 2), "f10"=matrix(c(-32, 32), nrow = 2),
                       "f11"=matrix(c(-600, 600), nrow = 2), "f12"=matrix(c(-50, 50), nrow = 2), 
                       "f13"=matrix(c(-50, 50), nrow = 2))
    iterasi <- input$iter[index]
    populasi <- input$popu[index]
    # test function
    tic("finished in")
    x <- metaAlgorithm(FUN = fungsi, optimType = optimType, numVar = numVar, numPopulation = populasi, 
                       maxIter = iterasi, rangeVar = rangeVar)
    time <- toc()
    print(x)
    time <- time$toc - time$tic
    optimum <- calcFitness(fungsi, 1, t(x))
    
    # save result
    x <- matrix(x, nrow = 1)
    fungsi <- input$fungsi[index]
    result <- data.frame(fungsi, iterasi, populasi, x, optimum, time)
    output <- rbind(output, result)
  }
  write.csv2(output, file = outputFile)
}

experiment(CLONALG, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultCLONALG.csv")
experiment(DE, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultDE.csv")
experiment(SFL, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultSFL.csv")
experiment(CSO, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultCSO.csv")
experiment(ABC, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultABC.csv")
experiment(KH, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultKH.csv")
experiment(CS, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultCS.csv")
set.seed(30)
experiment(BA, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultBA.csv")
experiment(GBS, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultGBS.csv")
experiment(BHO, inputFile = "./testFunction/test.csv", outputFile = "./testFunction/resultBHO.csv")
