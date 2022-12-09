# Create a transition matrix of a random walk on a square lattice graph
# This transition matrix will be used to calculate the order of the expected
# time to consensus using a formula from Cooper, et al, paper.
library(MASS)

# NOTE: the transition matrix returned is NOT SCALED
# to scale, add: tmat <- fractions(tmat/rowSums(tmat))
transitionMatrix.RandomWalk.lattice <- function(row, col){
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
tmat <- transitionMatrix.RandomWalk.lattice(3,2)
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

#-------#
# Looking at the complete graph

# random walk on a complete graph
# notice that the transition matrix is a matrix of 1's with 0's on the diagonal
transitionMatrix.RandomWalk.complete <- function(N){
  tmat <- matrix(data = 0, nrow = N, ncol = N) # initialize transition matrix
  
  for(i in 1:N){
    for(j in 1:N){
      if(j != i){
        tmat[i,j] <- 1
      }
    }
  }
  
  return(tmat)
}


tmat <- transitionMatrix.RandomWalk.complete(6)
N <- 6

## calculcate 2m
two.m <- sum(tmat)

## calculate nu
d2n <- two.m^2/N
nu <- sum(rowSums(tmat)^2) / d2n

## extract second largest eigenvalue
tmat <- fractions(tmat/rowSums(tmat))
eigen(tmat)$values # second largest eigenvalue: 0.2
lambda <- 0.2

## calculate order
N/(nu*(1-lambda)) # 7.5

## empirical result: 15.5
mean(VM_complete(N,10000))
