# prereqs
library(ggplot2)
library(gganimate)
library(gifski)
library(MASS)
library(magick)
theme_set(theme_bw())

# Uniformly selects neighbor of an individual using classic movements
# Parameters
  # nbhd: Matrix representing the different individuals
  # i: Row number of individual
  # j: Column number of individual
  # boundaries: If TRUE, we assume boundaries, otherwise edges and corners are
  #             linked to the other side of the matrix
  #             Default is set to TRUE
# Returns: row and column number of a randomly selected neighbor
neighbor <- function(nbhd,i,j,boundaries = TRUE){
  a <- 0
  b <- 0
  
  # possible movements (up/down/left/right)
  movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
  
  if(boundaries){
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
  } else {
    if(a == 0){
      a <- nrow(nbhd)
    } 
    if(a == nrow(nbhd) + 1){
      a <- 1
    }
    if(b == 0){
      b <- ncol(nbhd)
    }
    if(b == ncol(nbhd) + 1){
      b <- 1
    }
  }
  
  
  return(c(a,b))
}

# Classic voter model, used in VM theory.R file
# Assumes boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Time to consensus of each observation
VM <- function(row, col, s, pOne = 0.5){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne,pOne))
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(c(1:row), 1)
      j <- sample(c(1:col), 1)
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      a <- lst[1]
      b <- lst[2]
      
      nbhd[i,j] <- nbhd[a,b]
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}


VM_func <- function(fixed.row, col.range){
  avg <- numeric()
  var <- numeric()
  sec <- numeric()
  j <- 1
  for(i in col.range){
    print(i)
    time <- VM(fixed.row, i, 100)
    avg[j] <- mean(time)
    var[j] <- var(time)
    sec[j] <- mean(time^2)
    j <- j + 1
  }
  
  return(c(avg, var, sec))
}

# ----- #

# Confidence voter model - marginal
# Assumes boundaries
# Parameters:
  # row: number of rows
  # col: number of columns
  # s: number of observations
# Returns: Time to consensus of each observation
CVM_marg <- function(row, col, s){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    states <- sample(c(replicate(ceiling(N/2),'c'),replicate(ceiling(N/2), 'u')),N,replace = FALSE)
    confidence <- matrix(data = states, nrow = row, ncol = col)
      
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      neighb <- nbhd[lst[1],lst[2]]
      
      
      if(indiv == neighb){ # indiv becomes confident
        confidence[i,j] <- 'c'
      } else{ 
        if(confidence[i,j] == 'c'){ # indiv becomes unsure
          confidence[i,j] <- 'u'
        } else{ # indiv changes state
          nbhd[i,j] <- neighb
        }
      }
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Confidence voter model - extremal
# Assumes boundaries
# Parameters:
  # row: number of rows
  # col: number of columns
  # s: number of observations
# Returns: Time to consensus of each observation
CVM_extr <- function(row, col, s){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    states <- sample(c(replicate(ceiling(N/2),'c'),replicate(ceiling(N/2), 'u')),N,replace = FALSE)
    confidence <- matrix(data = states, nrow = row, ncol = col)
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      neighb <- nbhd[lst[1],lst[2]]
      
      
      if(indiv == neighb){ # indiv becomes confident
        confidence[i,j] <- 'c'
      } else{ 
        if(confidence[i,j] == 'c'){ # indiv becomes unsure
          confidence[i,j] <- 'u'
        } else{ # indiv changes state and becomes confidence in new state
          nbhd[i,j] <- neighb
          confidence[i,j] <- 'c'
        }
      }
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Models the time to consensus of the confidence voter models as a function 
# over a number of columns
# Parameters:
  # fixed.row: Number of rows that will remain fixed through each simulation
  # col.range: Range of columns. Each will be used twice: once in marginal, and 
  #            once in extremal case.
# Returns: Average time to, variance of, second moment of each dimension, for 
#          both marginal and extremal cases (6 vectors total)
CVM_func <- function(fixed.row, col.range){
  avg_m <- numeric()
  avg_e <- numeric()
  var_m <- numeric()
  var_e <- numeric()
  sec_m <- numeric()
  sec_e <- numeric()
  
  for(i in col.range){
    print(i)
    time_m <- CVM_marg(fixed.row, i, 100)
    time_e <- CVM_extr(fixed.row, i, 100)
    avg_m[i] <- mean(time_m)
    avg_e[i] <- mean(time_e)
    
    var_m[i] <- var(time_m)
    var_e[i] <- var(time_e)
    
    sec_m[i] <- mean(time_m^2)
    sec_e[i] <- mean(time_e^2)
  }
  
  return(c(avg_m, var_m, sec_m, avg_e, var_e, sec_e))
}

# ----- #

# Runs the same as VM, but
# Assume NO boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
# Returns: Time to consensus of each observation
VM_noBound <- function(row,col,s){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j, boundaries = F)
      a <- lst[1]
      b <- lst[2]
      
      nbhd[i,j] <- nbhd[a,b]
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Models the time to consensus of VM_noBound as a function over a number of 
# columns
# Parameters:
  # fixed.row: Number of rows that will remain fixed through each simulation
  # col.range: Range of columns
# Returns: Average time to, variance of, second moment of each dimension 
#          (3 vectors total)
VM_func_noBound <- function(fixed.row, col.range){
  avg <- numeric()
  var <- numeric()
  sec <- numeric()
  
  for(i in col.range){
    time <- VM_noBound(fixed.row, i, 100000)
    avg <- append(avg,mean(time))
    var <- append(var,var(time))
    sec <- append(sec,mean(time^2))
  }
  
  return(c(avg, var, sec))
}

# ----- #

# Confidence voter model - marginal (Version 2)
# Runs similarly to CVM_marg, but state/confidence only change if neighbor is 
# confident. 
# We can also adjust the probability that an individual is convinced by a 
# neighbor (conditioned on that neighbor being confident).
# Assumes boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
  # change: Probability that an individual will be convinced by a confident 
  #         neighbor (constant for each individual)
  #         Default is 1
# Returns: Time to consensus of each observation
CVM_marg_V2 <- function(row,col,s,change = 1){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    states <- sample(c(replicate(ceiling(N/2),'c'),replicate(ceiling(N/2), 'u')),N,replace = FALSE)
    confidence <- matrix(data = states, nrow = row, ncol = col)
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      neighb <- nbhd[lst[1],lst[2]]
      
      if(runif(1) <= change){
        if(indiv == neighb & confidence[lst[1],lst[2]] == 'c'){ # same opinion & neighbor is confidence
          confidence[i,j] <- 'c'
        } else{ 
          if(confidence[i,j] == 'c' & confidence[lst[1],lst[2]] == 'c'){ # indiv becomes unsure
            confidence[i,j] <- 'u'
          } else{ # indiv changes state
            if(indiv != neighb & confidence[i,j] == 'u'){
              if(confidence[lst[1],lst[2]] == 'c'){
                nbhd[i,j] <- neighb
              }
            }
          }
        }
      }
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Confidence voter model - extremal (Version 2)
# Runs similarly to CVM_extr, but state/confidence only change if neighbor is 
# confident. 
# We can also adjust the probability that an individual is convinced by a 
# neighbor (conditioned on that neighbor being confident).
# Assumes boundaries
# Parameters:
# row: Number of rows
# col: Number of columns
# s: Number of observations
# change: Probability that an individual will be convinced by a confident 
#         neighbor (constant for each individual)
#         Default is 1
# Returns: Time to consensus of each observation
CVM_extr_V2 <- function(row, col, s, change = 1){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    states <- sample(c(replicate(ceiling(N/2),'c'),replicate(ceiling(N/2), 'u')),N,replace = FALSE)
    confidence <- matrix(data = states, nrow = row, ncol = col)
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      neighb <- nbhd[lst[1],lst[2]]
      
      if(runif(1) <= change){
        if(indiv == neighb & confidence[lst[1],lst[2]] == 'c'){ # indiv becomes confident
          confidence[i,j] <- 'c'
        } else{ 
          if(confidence[i,j] == 'c' & confidence[lst[1],lst[2]] == 'c'){ # indiv becomes unsure
            if(indiv != neighb){
              confidence[i,j] <- 'u'
            }
          } else{ # indiv changes state and becomes confidence in new state
            if(confidence[i,j] == 'u' & confidence[lst[1],lst[2]] == 'c'){
              if(indiv != neighb){
                nbhd[i,j] <- neighb
                confidence[i,j] <- 'c'
              }
            }
          }
        }
      }
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Models the time to consensus of the confidence voter models (Version 2) as a  
# function over a number of columns
# Parameters:
  # fixed.row: Number of rows that will remain fixed through each simulation
  # col.range: Range of columns. Each will be used twice: once in marginal, and 
  #            once in extremal case.
  # change: Probability that an individual will be convinced by a confident 
  #         neighbor (constant for each individual)
  #         Default is 1
# Returns: Average time to, variance of, second moment of each dimension, for 
#          both marginal and extremal cases (6 vectors total)
CVM_func_V2 <- function(fixed.row, col.range, change = 1){
  avg_m <- numeric()
  avg_e <- numeric()
  var_m <- numeric()
  var_e <- numeric()
  sec_m <- numeric()
  sec_e <- numeric()
  
  for(i in col.range){
    print(i)
    time_m <- CVM_marg_V2(fixed.row, i, 1000, change = change)
    time_e <- CVM_extr_V2(fixed.row, i, 1000, change = change)
    avg_m[i] <- mean(time_m)
    avg_e[i] <- mean(time_e)
    
    var_m[i] <- var(time_m)
    var_e[i] <- var(time_e)
    
    sec_m[i] <- mean(time_m^2)
    sec_e[i] <- mean(time_e^2)
  }
  
  return(c(avg_m, var_m, sec_m, avg_e, var_e, sec_e))
}

# ----- #

# animating the classic voter model: special case 1
# for fun
VM_anim_special_1 <- function(row = 8, col = 12){
  
  N <- row*col
  k <- 1
  # initial state will not be random
  nbhd <- matrix(data = replicate(row*col, 0), nrow = row, ncol = col)
  nbhd[3:6,5:8] <- 1
  
  nbhd.df <- data.frame(reshape2::melt(nbhd), k = replicate(N,k))
  
  while(TRUE){
    k <- k + 1  
    if(sum(nbhd) == N | sum(nbhd) == 0){
      break;
    }
    
    # randomly select an individual
    i <- round(runif(1,1,row))
    j <- round(runif(1,1,col))
    
    indiv <- nbhd[i,j]
    
    # select a neighbor
    lst <- neighbor(nbhd, i, j)
    a <- lst[1]
    b <- lst[2]
    
    nbhd[i,j] <- nbhd[a,b]
    
    nbhd.df <- rbind(nbhd.df, data.frame(reshape2::melt(nbhd), k = replicate(N, k)))
  }
  
  anim <- ggplot(data = nbhd.df, mapping = aes(x = Var2, y = Var1)) + 
    geom_point(aes(size = 2, color = value)) +
    scale_x_continuous(name = "column number", breaks = 1:col) +
    scale_y_continuous(name = "row number", breaks = 1:row) +
    scale_colour_gradient(low = "yellow", high = "blue") +
    transition_time(nbhd.df$k)
  
  animate(anim, renderer = gifski_renderer("Classic_Special_Case_1.gif", loop = FALSE))
}

# animating the classic voter model: special case 2
# for fun
VM_anim_special_2 <- function(row = 8, col = 12){
  
  N <- row*col
  k <- 1
  # initial state will be approx. half random half manually added
  states <- round(runif(N))
  nbhd <- matrix(data = states, nrow = row, ncol = col)
  nbhd[1:4, 9:12] <- 0 # 5x5 in top right
  nbhd[5:8, 1:4]  <- 1 # 5x5 in bottom left
  
  nbhd.df <- data.frame(reshape2::melt(nbhd), k = replicate(N,k))
  
  while(TRUE){
    k <- k + 1  
    if(sum(nbhd) == N | sum(nbhd) == 0){
      break;
    }
    
    # randomly select an individual
    i <- round(runif(1,1,row))
    j <- round(runif(1,1,col))
    
    indiv <- nbhd[i,j]
    
    # select a neighbor
    lst <- neighbor(nbhd, i, j)
    a <- lst[1]
    b <- lst[2]
    
    nbhd[i,j] <- nbhd[a,b]
    
    nbhd.df <- rbind(nbhd.df, data.frame(reshape2::melt(nbhd), k = replicate(N, k)))
  }
  
  anim <- ggplot(data = nbhd.df, mapping = aes(x = Var2, y = Var1)) + 
    geom_point(aes(size = 2, color = value)) +
    scale_x_continuous(name = "column number", breaks = 1:col) +
    scale_y_continuous(name = "row number", breaks = 1:row) +
    scale_colour_gradient(low = "yellow", high = "blue") +
    transition_time(nbhd.df$k)
  
  animate(anim, renderer = gifski_renderer("Classic_Special_Case_2.gif", loop = FALSE))
}

# ----- #

# Runs the same as VM_noBound, so
# Assumes NO boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Vector of which value wins consensus in each simulation
VM_win <- function(row,col,s, pOne = 0.5){
  
  N <- row*col
  winner <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne, pOne))
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    
    while(TRUE){
      
      if( sum(nbhd) == N ){
        winner[k] <- 1
        break
      }
      
      if( sum(nbhd) == 0 ){
        winner[k] <- 0
        break
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j, boundaries = F)
      a <- lst[1]
      b <- lst[2]
      
      nbhd[i,j] <- nbhd[a,b]
    
    }
  }
  print(paste("One wins at probability ", pOne,": ", sum(winner == 1)/s, sep = ""))
  return(winner)
}

# Runs the same as VM_noBound, so
# Assumes NO boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Percent of individuals exist with state '1' at each time
VM_percents <- function(row,col, pOne = 0.5){
  
  N <- row*col
  percentOne <- numeric()
  states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne, pOne))
  nbhd <- matrix(data = states, nrow = row, ncol = col)
  percentOne[1] <- sum(nbhd == 1)/N
  
  while(TRUE){
    
    if( sum(nbhd) == N ){
      break
    }
    
    if( sum(nbhd) == 0 ){
      break
    }
    
    # randomly select an individual
    i <- sample(1:row, 1)
    j <- sample(1:col, 1)
    
    indiv <- nbhd[i,j]
    
    # select a neighbor
    lst <- neighbor(nbhd, i, j, boundaries = F)
    a <- lst[1]
    b <- lst[2]
    
    nbhd[i,j] <- nbhd[a,b]
    
    percentOne <- append(percentOne, sum(nbhd == 1)/N)
  }
  
  return(percentOne)
}

# ----- #

# Individuals may now interact with neighbors "diagonal" to them
# Uses an internal neighbor function, which includes both diagonal movements
# and classic movements
# Assumes NO boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Time to consensus of each observation
VM_diag <- function(row, col, s, pOne = 0.5){
  
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0),
                      c(1,1),c(1,-1),c(-1,1),c(-1,-1))
    
    
    index <- sample(1:8,1, replace = T)
    a <- i + movements[[index]][1]
    b <- j + movements[[index]][2]
    
    if(a == 0){
      a <- nrow(nbhd)
    } 
    if(a == nrow(nbhd) + 1){
      a <- 1
    }
    if(b == 0){
      b <- ncol(nbhd)
    }
    if(b == ncol(nbhd) + 1){
      b <- 1
    }
    
    return(c(a,b))
  }
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne, pOne))
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(c(1:row), 1)
      j <- sample(c(1:col), 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      a <- lst[1]
      b <- lst[2]
      
      nbhd[i,j] <- nbhd[a,b]
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Models the time to consensus of VM_diag as a function over a number of 
# columns
# Parameters:
  # fixed.row: Number of rows that will remain fixed through each simulation
  # col.range: Range of columns
# Returns: Average time to, variance of, second moment of each dimension 
#          (3 vectors total)
VM_diag_func <- function(fixed.row, col.range){
  avg <- numeric()
  var <- numeric()
  sec <- numeric()
  
  for(i in col.range){
    print(i)
    time <- VM_diag(fixed.row, i, 100000)
    avg <- append(avg,mean(time))
    var <- append(var,var(time))
    sec <- append(sec,mean(time^2))
  }
  
  return(c(avg, var, sec))
}

# ----- #

# Confidence voter model - marginal (Version 3)
# Each individual will have their own probability of being influenced by a 
# confident voter
# Assumes boundaries, and only classic movements
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Time to consensus of each observation
CVM_marg_V3 <- function(row, col, s, pOne = 0.5){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE, prob = c(1 - pOne, pOne))
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    states <- sample(c(replicate(ceiling(N/2),'c'),replicate(ceiling(N/2), 'u')),N,replace = FALSE)
    confidence <- matrix(data = states, nrow = row, ncol = col)
    change <- matrix(data = runif(N), nrow = row, ncol = col) # probability that an individual is willing to switch
                                                              # opinions after interacting with a confident neighbor
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      neighb <- nbhd[lst[1],lst[2]]
      
      if(runif(1) < change[i,j]){  
        if(indiv == neighb & confidence[lst[1],lst[2]] == 'c'){ # same opinion & neighbor is confidence
          confidence[i,j] <- 'c'
        } else{ 
          if(confidence[i,j] == 'c' & confidence[lst[1],lst[2]] == 'c'){ # indiv becomes unsure
            confidence[i,j] <- 'u'
          } else{ # indiv changes state
            if(indiv != neighb & confidence[i,j] == 'u'){
              if(confidence[lst[1],lst[2]] == 'c'){
                nbhd[i,j] <- neighb
              }
            }
          }
        }
      }
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Confidence voter model - extremal (Version 3)
# Each individual will have their own probability of being influenced by a 
# confident voter
# Assumes boundaries, and only classic movements
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Time to consensus of each observation
CVM_extr_V3 <- function(row, col, s, pOne = 0.5){
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    states <- sample(c(replicate(ceiling(N/2),'c'),replicate(ceiling(N/2), 'u')),N,replace = FALSE)
    confidence <- matrix(data = states, nrow = row, ncol = col)
    change <- matrix(data = runif(N), nrow = row, ncol = col) # probability that an individual is willing to switch
                                                              # opinions after interacting with a confident neighbor
    
    while(TRUE){
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break;
      }
      
      # randomly select an individual
      i <- sample(1:row, 1)
      j <- sample(1:col, 1)
      
      indiv <- nbhd[i,j]
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      neighb <- nbhd[lst[1],lst[2]]
      
      if(runif(1) < change[i,j]){
        if(indiv == neighb & confidence[lst[1],lst[2]] == 'c'){ # indiv becomes confident
          confidence[i,j] <- 'c'
        } else{ 
          if(confidence[i,j] == 'c' & confidence[lst[1],lst[2]] == 'c'){ # indiv becomes unsure
            if(indiv != neighb){
              confidence[i,j] <- 'u'
            }
          } else{ # indiv changes state and becomes confidence in new state
            if(confidence[i,j] == 'u' & confidence[lst[1],lst[2]] == 'c'){
              if(indiv != neighb){
                nbhd[i,j] <- neighb
                confidence[i,j] <- 'c'
              }
            }
          }
        }
      }
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# Models the time to consensus of the confidence voter models (Version 3) as a  
# function over a number of columns
# Parameters:
  # fixed.row: Number of rows that will remain fixed through each simulation
  # col.range: Range of columns. Each will be used twice: once in marginal, and 
  #            once in extremal case.
# Returns: Average time to, variance of, second moment of each dimension, for 
#          both marginal and extremal cases (6 vectors total)
CVM_func_V3 <- function(fixed.row, col.range){
  avg_m <- numeric()
  avg_e <- numeric()
  var_m <- numeric()
  var_e <- numeric()
  sec_m <- numeric()
  sec_e <- numeric()
  
  for(i in col.range){
    print(i)
    time_m <- CVM_marg_V3(fixed.row, i, 1000)
    time_e <- CVM_extr_V3(fixed.row, i, 1000)
    avg_m[i] <- mean(time_m)
    avg_e[i] <- mean(time_e)
    
    var_m[i] <- var(time_m)
    var_e[i] <- var(time_e)
    
    sec_m[i] <- mean(time_m^2)
    sec_e[i] <- mean(time_e^2)
  }
  
  
  
  return(c(avg_m, var_m, sec_m, avg_e, var_e, sec_e))
}

# uses internal neighbor methods (no boundaries)
# for fun
couplingVM <- function(row, col, pOne = 0.5){
  
  # select a neighbor using classic method: up/down/left/right
  neighborC <- function(nbhd,i,j){
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    index <- sample(1:4,1, replace = T)
    
    # coordinates of neighbor
    a <- i + movements[[index]][1]
    b <- j + movements[[index]][2]
    
    # remove boundaries
    if(a == 0){
      a <- nrow(nbhd)
    } 
    if(a == nrow(nbhd) + 1){
      a <- 1
    }
    if(b == 0){
      b <- ncol(nbhd)
    }
    if(b == ncol(nbhd) + 1){
      b <- 1
    }
    
    return(c(a,b))
  }
  
  # select a neighbor, but now, in addition to classic moves, we can move diagonally
  neighborD <- function(nbhd,i,j){
    
    # possible movements
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0),
                      c(1,1),c(1,-1),c(-1,1),c(-1,-1))
    index <- sample(1:8,1, replace = T)
    
    # coordinates of neighbor
    a <- i + movements[[index]][1]
    b <- j + movements[[index]][2]
    
    # remove boundaries
    if(a == 0){
      a <- nrow(nbhd)
    } 
    if(a == nrow(nbhd) + 1){
      a <- 1
    }
    if(b == 0){
      b <- ncol(nbhd)
    }
    if(b == ncol(nbhd) + 1){
      b <- 1
    }
    
    return(c(a,b))
  }
  
  # k <- 1 # for animation
  N <- row*col # total number of particles
  
  # percent of "1"s in each model
  percentOne.classic <- numeric()
  percentOne.diag <- numeric()
  
  # initialize each model with different set of states
  states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne, pOne))
  nbhd.classic <- matrix(data = states, nrow = row, ncol = col)
  
  states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne, pOne))
  nbhd.diag <- matrix(data = states, nrow = row, ncol = col)
  
  # record initial number of "1"s in each model
  percentOne.classic[1] <- sum(nbhd.classic == 1)/N
  percentOne.diag[1] <- sum(nbhd.diag == 1)/N
  
  # # record initial states for animation
  # nbhd.df.classic <- data.frame(reshape2::melt(nbhd.classic), k = replicate(N,k))
  # nbhd.df.diag <- data.frame(reshape2::melt(nbhd.classic), k = replicate(N,k))
  
  # run until either of the models reach consensus
  while(TRUE){
    k <- k + 1
    
    if(sum(nbhd.classic) == N | sum(nbhd.classic) == 0){
      break;
    }
    if(sum(nbhd.diag) == N | sum(nbhd.diag) == 0){
      break;
    }
    
    # randomly select an individual
    i <- sample(c(1:row), 1)
    j <- sample(c(1:col), 1)
    
    
    # select a neighbor
    lst <- neighborC(nbhd.classic, i, j)
    a <- lst[1]
    b <- lst[2]
    
    nbhd.classic[i,j] <- nbhd.classic[a,b]
    
    lst <- neighborD(nbhd.diag, i, j)
    a <- lst[1]
    b <- lst[2]
    
    nbhd.diag[i,j] <- nbhd.diag[a,b]
    
    percentOne.classic <- append(percentOne.classic, sum(nbhd.classic == 1)/N)
    percentOne.diag <- append(percentOne.diag, sum(nbhd.diag == 1)/N)
    
    # nbhd.df.classic <- rbind(nbhd.df.classic, data.frame(reshape2::melt(nbhd.classic), k = replicate(N, k)))
    # nbhd.df.diag <- rbind(nbhd.df.diag, data.frame(reshape2::melt(nbhd.diag), k = replicate(N, k)))
  }
  
  # # create side-by-side animation
  # anim1 <- ggplot(data = nbhd.df.classic, mapping = aes(x = Var2, y = Var1)) + 
  #   geom_point(aes(size = 2, color = value)) +
  #   scale_x_continuous(name = "column number", breaks = 1:col) +
  #   scale_y_continuous(name = "row number", breaks = 1:row) +
  #   scale_colour_gradient(low = "yellow", high = "blue") +
  #   labs(title = "Classic Movements") +
  #   transition_time(nbhd.df.classic$k)
  #  
  # anim2 <- ggplot(data = nbhd.df.diag, mapping = aes(x = Var2, y = Var1)) + 
  #   geom_point(aes(size = 2, color = value)) +
  #   scale_x_continuous(name = "column number", breaks = 1:col) +
  #   scale_y_continuous(name = "row number", breaks = 1:row) +
  #   scale_colour_gradient(low = "yellow", high = "blue") +
  #   labs(title = "Diagonal Movements") +
  #   transition_time(nbhd.df.diag$k)
  # 
  # anim.classic <- animate(anim1, nframes = 100, renderer = magick_renderer())
  # 
  # anim.diag <- animate(anim2, nframes = 100, renderer = magick_renderer())
  # 
  # final.anim <- image_append(c(anim.classic[1],anim.diag[1]))
  # 
  # for(i in 2:100){
  #   combined <- image_append(c(anim.classic[i],anim.diag[i]))
  #   final.anim <- append(final.anim,combined)
  # }
  #
  # anim_save("coupling.gif", final.anim)
  
  return(c(percentOne.classic, percentOne.diag))
}

# Modified Voter model - Version 4, used in VM theory.R file
# Here, the probability of being convinced is dependent on the proportion of
# agreeing voters. That is, if we have 3 individuals labeled "0" and 7 voters
# labeled "1", a "0" voter has a probability of 3/10 of being convinced, and a
# "1" voter has a probability of 7/10 of being convinced.
# We only run 1 observation
# Assumes boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Percent of population labeled "1" throughout the observation
V4_percent <- function(row, col){
  
  N <- row*col
  percentOne<- numeric()
  consensus <- FALSE
  
  states <- sample(c(0,1), N, replace = TRUE, prob = replicate(n,1/n))
  nbhd <- matrix(data = states, nrow = row, ncol = col)
  percentOne[1] <- sum(nbhd == 1)/N
  
  while(TRUE){
    pChange <- c(sum(nbhd == 0)/N,sum(nbhd == 1)/N)
    
    if(sum(nbhd) == N | sum(nbhd) == 0){
      break
    }
    
    # randomly select an individual
    i <- sample(c(1:row), 1)
    j <- sample(c(1:col), 1)
    
    # select a neighbor
    lst <- neighbor(nbhd, i, j, boundaries = T)
    a <- lst[1]
    b <- lst[2]
    
    if(pChange[nbhd[i,j] + 1] > runif(1)){
      nbhd[i,j] <- nbhd[a,b]
    }
    
    percentOne <- append(percentOne, sum(nbhd == 1)/N)
  }
  
  
  return(percentOne)
}


# Modified Voter model - Version 4, used in VM theory.R file
# Here, the probability of being convinced is dependent on the proportion of
# agreeing voters. That is, if we have 3 individuals labeled "0" and 7 voters
# labeled "1", a "0" voter has a probability of 3/10 of being convinced, and a
# "1" voter has a probability of 7/10 of being convinced.
# Assumes boundaries
# Parameters:
  # row: Number of rows
  # col: Number of columns
  # s: Number of observations
  # pOne: Probability that a spot in the matrix is initialized with '1'
  #       Default is 0.5
# Returns: Time to consensus of each observation
V4_time <- function(row, col, s, pOne = 0.5){
  N <- row*col
  consensusT <- numeric(s)
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne,pOne))
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    
    while(TRUE){
      pChange <- c(sum(nbhd == 0)/N, sum(nbhd == 1)/N)
      
      if(sum(nbhd) == N | sum(nbhd) == 0){
        break
      }
      
      # randomly select an individual
      i <- sample(c(1:row), 1)
      j <- sample(c(1:col), 1)
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      a <- lst[1]
      b <- lst[2]
      
      if(pChange[nbhd[i,j] + 1] > runif(1)){
        nbhd[i,j] <- nbhd[a,b]
      }
      
      consensusT[k] <- consensusT[k] + 1 
    }
  }
  return(consensusT)
}

# Models the time to consensus of VM_4 as a function over a number of columns
# Parameters:
  # fixed.row: Number of rows that will remain fixed through each simulation
  # col.range: Range of columns
# Returns: Average time to, variance of, second moment of each dimension 
  #          (3 vectors total)
V4_func <- function(fixed.row, col.range){
  avg <- numeric()
  var <- numeric()
  sec <- numeric()
  
  for(i in col.range){
    print(i)
    time <- V4_time(fixed.row, i, 1000)
    avg[i] <- mean(time)
    var[i] <- var(time)
    sec[i] <- mean(time^2)
  }
  
  return(c(avg, var, sec))
}

oneWins <- function(row, col){
  N <- row*col
  condition <- FALSE
  while(!condition){
    states <- sample(c(replicate(N-1, 0),1), N, replace = FALSE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    consensusT <- 0
    percentOne <- numeric()
    percentOne[1] <- (sum(nbhd == 1)/N)
    k <- 2
    while(TRUE){
      
      if(sum(nbhd) == 0){
        break
      }
      if(sum(nbhd) == N){
        condition <- TRUE
        break
      }
      
      # randomly select an individual
      i <- sample(c(1:row), 1)
      j <- sample(c(1:col), 1)
      
      # select a neighbor
      lst <- neighbor(nbhd, i, j)
      a <- lst[1]
      b <- lst[2]
      
      nbhd[i,j] <- nbhd[a,b]
      
      percentOne[k] <- (sum(nbhd == 1)/N)
      consensusT <- consensusT + 1
      k <- k + 1
    }
  }
  
  return(percentOne)
}
