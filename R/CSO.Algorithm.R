# Cat Swarm Optimization (CSO)

source('./R/metaheuristic.FunctionCollection.R')

ratio <- 0.3
numPopulation <- 13

numSeeking <- ratio * numPopulation
numTracing <- (1 - ratio) * numPopulation

seeking <- rep("seeking", round(numSeeking))
tracing <- rep("tracing", round(numTracing))

sample(c(seeking, tracing), replace = FALSE)
