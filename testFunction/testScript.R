# Test for CLONALG------------------------------------------------------------------------
library("tictoc", lib.loc="~/R/win-library/3.5")
source("./testFunction/13TestFunction.R")
source("./R/CLONALG.Algorithm.R")

#f1
tic(msg = "Finished in")
CLONALG(FUN=f1, optimType=test.spec$optimType[1] , numVar=test.spec$numVar[1], 
        rangeVar=matrix(c(test.spec$lowerBound[1], test.spec$upperBound[1]), nrow = 2))
toc()
#f2
CLONALG(FUN=f2, optimType=test.spec$optimType[2] , numVar=test.spec$numVar[2], 
        rangeVar=matrix(c(test.spec$lowerBound[2], test.spec$upperBound[2]), nrow = 2))
#f3
CLONALG(FUN=f3, optimType=test.spec$optimType[3] , numVar=test.spec$numVar[3], 
        rangeVar=matrix(c(test.spec$lowerBound[3], test.spec$upperBound[3]), nrow = 2))
#f4
CLONALG(FUN=f4, optimType=test.spec$optimType[4] , numVar=test.spec$numVar[4], 
        rangeVar=matrix(c(test.spec$lowerBound[4], test.spec$upperBound[4]), nrow = 2))
#f5
CLONALG(FUN=f5, optimType=test.spec$optimType[5] , numVar=test.spec$numVar[5], 
        rangeVar=matrix(c(test.spec$lowerBound[5], test.spec$upperBound[5]), nrow = 2))
#f6
CLONALG(FUN=f6, optimType=test.spec$optimType[6] , numVar=test.spec$numVar[6], 
        rangeVar=matrix(c(test.spec$lowerBound[6], test.spec$upperBound[6]), nrow = 2))
#f7
CLONALG(FUN=f7, optimType=test.spec$optimType[7] , numVar=test.spec$numVar[7], 
        rangeVar=matrix(c(test.spec$lowerBound[7], test.spec$upperBound[7]), nrow = 2))
#f8
CLONALG(FUN=f8, optimType=test.spec$optimType[8] , numVar=test.spec$numVar[8], 
        rangeVar=matrix(c(test.spec$lowerBound[8], test.spec$upperBound[8]), nrow = 2))
#f9
CLONALG(FUN=f9, optimType=test.spec$optimType[9] , numVar=test.spec$numVar[9], 
        rangeVar=matrix(c(test.spec$lowerBound[9], test.spec$upperBound[9]), nrow = 2))
#f10
CLONALG(FUN=f10, optimType=test.spec$optimType[10] , numVar=test.spec$numVar[10], 
        rangeVar=matrix(c(test.spec$lowerBound[10], test.spec$upperBound[10]), nrow = 2))
#f11
CLONALG(FUN=f11, optimType=test.spec$optimType[11] , numVar=test.spec$numVar[11], 
        rangeVar=matrix(c(test.spec$lowerBound[11], test.spec$upperBound[11]), nrow = 2))
#f12
CLONALG(FUN=f12, optimType=test.spec$optimType[12] , numVar=test.spec$numVar[12], 
        rangeVar=matrix(c(test.spec$lowerBound[12], test.spec$upperBound[12]), nrow = 2))
#f13
CLONALG(FUN=f13, optimType=test.spec$optimType[13] , numVar=test.spec$numVar[13], 
        rangeVar=matrix(c(test.spec$lowerBound[13], test.spec$upperBound[13]), nrow = 2))

# Test for BHO -----------------------------------------------------------------------------
source("./testFunction/13TestFunction.R")
source("./R/BHO.Algorithm.R")

#f1
BHO(FUN=f1, optimType=test.spec$optimType[1] , numVar=test.spec$numVar[1], 
        rangeVar=matrix(c(test.spec$lowerBound[1], test.spec$upperBound[1]), nrow = 2))
#f2
BHO(FUN=f2, optimType=test.spec$optimType[2] , numVar=test.spec$numVar[2], 
        rangeVar=matrix(c(test.spec$lowerBound[2], test.spec$upperBound[2]), nrow = 2))
#f3
BHO(FUN=f3, optimType=test.spec$optimType[3] , numVar=test.spec$numVar[3], 
        rangeVar=matrix(c(test.spec$lowerBound[3], test.spec$upperBound[3]), nrow = 2))
#f4
BHO(FUN=f4, optimType=test.spec$optimType[4] , numVar=test.spec$numVar[4], 
        rangeVar=matrix(c(test.spec$lowerBound[4], test.spec$upperBound[4]), nrow = 2))
#f5
BHO(FUN=f5, optimType=test.spec$optimType[5] , numVar=test.spec$numVar[5], 
        rangeVar=matrix(c(test.spec$lowerBound[5], test.spec$upperBound[5]), nrow = 2))
#f6
BHO(FUN=f6, optimType=test.spec$optimType[6] , numVar=test.spec$numVar[6], 
        rangeVar=matrix(c(test.spec$lowerBound[6], test.spec$upperBound[6]), nrow = 2))
#f7
BHO(FUN=f7, optimType=test.spec$optimType[7] , numVar=test.spec$numVar[7], 
        rangeVar=matrix(c(test.spec$lowerBound[7], test.spec$upperBound[7]), nrow = 2))
#f8
BHO(FUN=f8, optimType=test.spec$optimType[8] , numVar=test.spec$numVar[8], 
        rangeVar=matrix(c(test.spec$lowerBound[8], test.spec$upperBound[8]), nrow = 2))
#f9
BHO(FUN=f9, optimType=test.spec$optimType[9] , numVar=test.spec$numVar[9], 
        rangeVar=matrix(c(test.spec$lowerBound[9], test.spec$upperBound[9]), nrow = 2))
#f10
BHO(FUN=f10, optimType=test.spec$optimType[10] , numVar=test.spec$numVar[10], 
        rangeVar=matrix(c(test.spec$lowerBound[10], test.spec$upperBound[10]), nrow = 2))
#f11
BHO(FUN=f11, optimType=test.spec$optimType[11] , numVar=test.spec$numVar[11], 
        rangeVar=matrix(c(test.spec$lowerBound[11], test.spec$upperBound[11]), nrow = 2))
#f12
BHO(FUN=f12, optimType=test.spec$optimType[12] , numVar=test.spec$numVar[12], 
        rangeVar=matrix(c(test.spec$lowerBound[12], test.spec$upperBound[12]), nrow = 2))
#f13
BHO(FUN=f13, optimType=test.spec$optimType[13] , numVar=test.spec$numVar[13], 
        rangeVar=matrix(c(test.spec$lowerBound[13], test.spec$upperBound[13]), nrow = 2))

# Test for DE -----------------------------------------------------------------------------
source("./testFunction/13TestFunction.R")
source("./R/DE.Algorithm.R")

#f1
DE(FUN=f1, optimType=test.spec$optimType[1] , numVar=test.spec$numVar[1], 
    rangeVar=matrix(c(test.spec$lowerBound[1], test.spec$upperBound[1]), nrow = 2))
#f2
DE(FUN=f2, optimType=test.spec$optimType[2] , numVar=test.spec$numVar[2], 
    rangeVar=matrix(c(test.spec$lowerBound[2], test.spec$upperBound[2]), nrow = 2))
#f3
DE(FUN=f3, optimType=test.spec$optimType[3] , numVar=test.spec$numVar[3], 
    rangeVar=matrix(c(test.spec$lowerBound[3], test.spec$upperBound[3]), nrow = 2))
#f4
DE(FUN=f4, optimType=test.spec$optimType[4] , numVar=test.spec$numVar[4], 
    rangeVar=matrix(c(test.spec$lowerBound[4], test.spec$upperBound[4]), nrow = 2))
#f5
DE(FUN=f5, optimType=test.spec$optimType[5] , numVar=test.spec$numVar[5], 
    rangeVar=matrix(c(test.spec$lowerBound[5], test.spec$upperBound[5]), nrow = 2))
#f6
DE(FUN=f6, optimType=test.spec$optimType[6] , numVar=test.spec$numVar[6], 
    rangeVar=matrix(c(test.spec$lowerBound[6], test.spec$upperBound[6]), nrow = 2))
#f7
DE(FUN=f7, optimType=test.spec$optimType[7] , numVar=test.spec$numVar[7], 
    rangeVar=matrix(c(test.spec$lowerBound[7], test.spec$upperBound[7]), nrow = 2))
#f8
DE(FUN=f8, optimType=test.spec$optimType[8] , numVar=test.spec$numVar[8], 
    rangeVar=matrix(c(test.spec$lowerBound[8], test.spec$upperBound[8]), nrow = 2))
#f9
DE(FUN=f9, optimType=test.spec$optimType[9] , numVar=test.spec$numVar[9], 
    rangeVar=matrix(c(test.spec$lowerBound[9], test.spec$upperBound[9]), nrow = 2))
#f10
DE(FUN=f10, optimType=test.spec$optimType[10] , numVar=test.spec$numVar[10], 
    rangeVar=matrix(c(test.spec$lowerBound[10], test.spec$upperBound[10]), nrow = 2))
#f11
DE(FUN=f11, optimType=test.spec$optimType[11] , numVar=test.spec$numVar[11], 
    rangeVar=matrix(c(test.spec$lowerBound[11], test.spec$upperBound[11]), nrow = 2))
#f12
DE(FUN=f12, optimType=test.spec$optimType[12] , numVar=test.spec$numVar[12], 
    rangeVar=matrix(c(test.spec$lowerBound[12], test.spec$upperBound[12]), nrow = 2))
#f13
DE(FUN=f13, optimType=test.spec$optimType[13] , numVar=test.spec$numVar[13], 
    rangeVar=matrix(c(test.spec$lowerBound[13], test.spec$upperBound[13]), nrow = 2))
