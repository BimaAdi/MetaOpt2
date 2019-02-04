# 13 Test Function
test.spec <- data.frame(id = c("f1","f2","f3","f4","f5","f6","f7","f8","f9","f10","f11","f12","f13"),
                            name = c("Sphere Model", "Schewefel's Problem 2.22", "Schewefel's Problem 1.2",
                                      "Schewefel's Problem 2.21", "Generalized Rosenbrock's", "Step Function",
                                      "Quartic with Noise", "Generalized Schewefel's Problem 2.26", "Generalized Rastrigin's Function",
                                      "Ackley's Function", "Generalized Griewank Function", "Generalized Penalized Function 1",
                                      "Generalized Penalized Function 2"),
                            optimType=rep("MIN", 13),
                            numVar=c(1, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
                            lowerBound = c(-100, -10, -100, -100, -30, -100, -1.28, -500, -5.12, -32, -600, -50, -50),
                            upperBound = c(100, 10, 100, 100, 30, 100, 1.28, 500, 5.12, 32, 600, 50, 50),
                            x = c(0, 0, 0, 0, 1, 0, 0, 420.9687, 0, 0, 0, 1, 1),
                            fmin = c(0, 0, 0, 0, 0, 0, 0, -Inf, 0, 0, 0, 0, 0))
f1 <- function(X){
  return(sum(X^2))
}
f2 <- function(x){
  return(sum(abs(x)+prod(abs(x))))
}
f3 <- function(x){
  dim <- length(x)
  result <- 0
  for(i in 1:dim){
    result <- result + sum(x[1:i])^2
  }
  return(result)
}
f4 <- function(x){
  return(max(abs(x)))
}
f5 <- function(x){
  dim <- length(x)
  result <- sum(100*(x[2:dim]-(x[1:dim-1]^2))^2+(x[1:dim-1]-1)^2)
  return(result)
}
f6 <- function(x){
  result <- sum(abs((x+0.5))^2)
  return(result)
}
f7 <- function(x){
  dim <- length(x)
  result <- sum(c(1:dim)*(x^4))+runif(1)
  return(result)
}
f8 <- function(x){
  result <- sum(-x*sin(sqrt(abs(x))))
  return(result)
}
f9 <- function(x){
  dim <- length(x)
  result <- sum(x^2-10*cos(2*pi*x))+10*dim
  return(result)
}
f10 <- function(x){
  dim <- length(x)
  result <- -20*exp(-0.2*sqrt(sum(x^2)/dim))-exp(sum(cos(2*pi*x))/dim)+20+exp(1)
  return(result)
}
f11 <- function(x){
  dim <- length(x)
  result <- sum(x^2)/4000-prod(cos(x/sqrt(c(1:dim))))+1
  return(result)
}
Ufun <- function(x,a,k,m){
  result <- k*((x-a)^m)*(x>a)+k*((-x-a)^m)*(x<(-a))
  return(result)
}
f12 <- function(x){
  dim <- length(x)
  result <- (pi/dim)*(10*((sin(pi*(1+(x[1]+1)/4)))^2)+sum((((x[1:dim-1]+1)/4)^2)*(1+10*((sin(pi*(1+(x[2:dim]+1)/4))))^2))+((x[dim]+1)/4)^2)+sum(Ufun(x,10,100,4))
  return(result)
}
f13 <- function(x){
  dim <- length(x)
  result <- 0.1*((sin(3*pi*x[1]))^2+sum((x[1:dim-1]-1)^2*(1+(sin(3*pi*x[2:dim]))^2))+((x[dim]-1)^2)*(1+(sin(2*pi*x[dim]))^2))+sum(Ufun(x,5,100,4))
  return(result)
}