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