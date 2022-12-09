# This file will allow us to generate a transition matrix of any size model without
# calculating anything by hand.

# Prerequisites
source("Voter-model-simulations.R")

# Find all possible states of a nxp model
# Parameters:
  # col: Number of columns
  # row: Number of rows
# Returns a matrix, where each row is a possible state that the model can take
states <- function(row, col){
  N <- col*row
  possibilities <- matrix(data = c(rep(0,N),rep(1,N)), ncol = N, nrow = 2, byrow = TRUE)
  
  while(nrow(possibilities) < 2^N){
    if(nrow(possibilities) %% 1000 == 0){
      print(nrow(possibilities))
    }
    temp <- sample(c(0,1), N, replace = TRUE)

    # check if state already exists in data frame
    
    tryCatch(
      {
        # check if state already exists
        if(which(rowSums(t(t(possibilities) == temp)) == N) > 0){
          next
        }
      },
      error = function(cond){ # logical error: new state found
        # possibilities <- rbind(possibilities, temp)
      }
    )
    possibilities <- rbind(possibilities, temp)
    
  }
  
  row.names(possibilities) <- 1:(2^N)
  
  return(possibilities)
}

# Automatically calculates the transition matrix of a classic voter model
  # Parameters:
  # col: Number of columns
  # rows: Number of rows
# Returns: Transition matrix
transitionMatrix.Classic <- function(row,col){
  states.matrix <- states(row,col)
  nr <- nrow(states.matrix)
  N <- row*col
  n <- 1000000
  
  # find transition probabilities for first state to initialize t.matrix
  s0 <- matrix(data = states.matrix[1,], ncol = col, nrow = row, byrow = T)
  total <- numeric(nr)
  for(k in 1:n){
    temp <- s0
    # randomly select an individual
    i <- sample(c(1:row), 1)
    j <- sample(c(1:col), 1)
    
    # select a neighbor
    lst <- neighbor(s0, i, j)
    a <- lst[1]
    b <- lst[2]
    
    temp[i,j] <- s0[a,b]
    
    temp <- as.vector(temp)
    
    total[which(rowSums(t(t(states.matrix) == temp)) == N)] <- total[which(rowSums(t(t(states.matrix) == temp)) == N)] + 1
  }
  
  t.matrix <- matrix(data <- total/n, ncol = nr, nrow = 1, byrow = T)
  
  # find the rest of the probabilities
  for(l in 2:nr){
    print(l)
    s0 <- matrix(data = states.matrix[l,], ncol = col, nrow = row, byrow = T)
    total <- numeric(nr)
    for(k in 1:n){
      temp <- s0
      # randomly select an individual
      i <- sample(c(1:row), 1)
      j <- sample(c(1:col), 1)
      
      # select a neighbor
      lst <- neighbor(s0, i, j)
      a <- lst[1]
      b <- lst[2]
      
      temp[i,j] <- s0[a,b]
      
      temp <- as.vector(temp)
      
      total[which(rowSums(t(t(states.matrix) == temp)) == N)] <- total[which(rowSums(t(t(states.matrix) == temp)) == N)] + 1
    }
    t.matrix <- rbind(t.matrix, total/n)
  }
  
  return(t.matrix)
}

# Automatically calculates the transition matrix of a modified voter model
# Parameters:
  # col: Number of columns
  # rows: Number of rows
# Returns: Transition matrix
transitionMatrix.Modified <- function(row,col){
  states.matrix <- states(row,col)
  nr <- nrow(states.matrix)
  N <- row*col
  n <- 10000
  
  # find transition probabilities for first state to initialize t.matrix
  s0 <- matrix(data = states.matrix[1,], ncol = col, nrow = row, byrow = T)
  total <- numeric(nr)
  for(k in 1:n){
    temp <- s0
    pChange <- c(sum(temp == 0)/N, sum(temp == 1)/N)
    
    # randomly select an individual
    i <- sample(c(1:row), 1)
    j <- sample(c(1:col), 1)
    
    # select a neighbor
    lst <- neighbor(s0, i, j)
    a <- lst[1]
    b <- lst[2]
    
    if(pChange[s0[i,j] + 1] > runif(1)){
      temp[i,j] <- s0[a,b]
    }
    
    temp <- as.vector(temp)
    
    total[which(rowSums(t(t(states.matrix) == temp)) == N)] <- total[which(rowSums(t(t(states.matrix) == temp)) == N)] + 1
  }
  
  t.matrix <- matrix(data <- total/n, ncol = nr, nrow = 1, byrow = T)
  
  # find the rest of the probabilities
  for(l in 2:nr){
    print(l)
    s0 <- matrix(data = states.matrix[l,], ncol = col, nrow = row, byrow = T)
    total <- numeric(nr)
    for(k in 1:n){
      temp <- s0
      pChange <- c(sum(temp == 0)/N, sum(temp == 1)/N)
      
      # randomly select an individual
      i <- sample(c(1:row), 1)
      j <- sample(c(1:col), 1)
      
      # select a neighbor
      lst <- neighbor(s0, i, j)
      a <- lst[1]
      b <- lst[2]
      
      if(pChange[s0[i,j] + 1] > runif(1)){
        temp[i,j] <- s0[a,b]
      }
      
      temp <- as.vector(temp)
      
      total[which(rowSums(t(t(states.matrix) == temp)) == N)] <- total[which(rowSums(t(t(states.matrix) == temp)) == N)] + 1
    }
    t.matrix <- rbind(t.matrix, total/n)
  }
  
  return(t.matrix)
}

# Calculates the expected time to absorption
  # Parameters:
  # row: Number of rows
  # col: Number of columns
  # t.matrix: transition matrix
# Returns: Expected time to absorption
calculate <- function(t.matrix){
  nr <- nrow(t.matrix)
  nc <- ncol(t.matrix)
  Q <- t.matrix[3:nr,3:nc] # exclude absorbing states
  prob <- rep(1/(nr), nr-2)
  
  # expected time
  ET <- prob %*% (solve(diag(nr-2) - Q) %*% rep(1,nr-2))
  return(ET)
}

#  -- testing transitionMatrix.classic -- #
# calculate(transitionMatrix.Classic(2,2)) # 5.77 (empirical result is close to 5.9)

# -- for 3x2 matrix -- #
# calculate(transitionMatrix.Classic(3,2)) # 15.7 (by-hand computation: 18.31765)

#  -- testing transitionMatrix.Modified -- #
# calculate(transitionMatrix.Modified(2,2)) #  25.73097 (by-hand computation: 25.75)

# -- for 3x2 matrix -- #
# calculate(transitionMatrix.Modified(3,2)) # 134 (by-hand computation: 160.3006)