# Create a transition matrix of a random walk on a square lattice graph
# This transition matrix will be used to calculate the order of the expected
# time to consensus using a formula from [Author].
library(MASS)

# NOTE: the transition matrix returned is NOT SCALED
# to scale, add: tmat <- fractions(tmat/rowSums(tmat))
transitionMatrix.RandomWalk <- function(row, col){
  N <- row*col
  tmat <- matrix(data = 0, nrow = N, ncol = N) # initialize transition matrix
  RW <- matrix(data = 1:N, nrow = row, ncol = col,byrow = T) # initialize random walk
  
  for(i in 1:row){
    for(j in 1:col){
      if(j < col){
        tmat[RW[i,j],RW[i,j+1]] <- 1
      }
      if(j > 1){
        tmat[RW[i,j],RW[i,j-1]] <- 1
      }
      if(i > 1){
        tmat[RW[i,j],RW[i-1,j]] <- 1
      }
      if(i < row){
        tmat[RW[i,j],RW[i+1,j]] <- 1
      }
    }
  }
  
  return(tmat)
}

# Example calculation using cooper's equation: 3x2 model
tmat <- transitionMatrix.RandomWalk(3,2)
N <- 3*2

## calculcate 2m
two.m <- sum(tmat)

## calculate nu
d2n <- two.m^2/N
nu <- sum(rowSums(tmat)^2) / d2n

## extract second largest eigenvalue
tmat <- fractions(tmat/rowSums(tmat))
eigen(tmat)$values # second largest eigenvalue: 0.5
lambda <- 0.5

## calculate order
N/(nu*(1-lambda)) # 11.529

## empirical result: 18.1
