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
probMatrix <- apply(randomMatrix, c(1, 2), function(x){
  runif(1)
})
x[probMatrix < 0.5]   
col(x)[probMatrix < 0.5]

lowerBound[col(x)[probMatrix < 0.5]]
upperBound[col(x)[probMatrix < 0.5]]
