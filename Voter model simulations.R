## For each simulation, assume boundaries unless otherwise specifieds

# prereqs
library(ggplot2)
library(gganimate)
library(gifski)
theme_set(theme_bw())

# classic voter model, used in VM theory.R file
# row: number of rows
# col: number of columns
# s: number of observations
# pOne: probability that a spot in the matrix is initialized with '1'
VM <- function(row, col, s, pOne = 0.5){
  
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0

    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
        index <- sample(1:4,1, replace = T)
        a <- i + movements[[index]][1]
        b <- j + movements[[index]][2]
        
        if( a > nrow(nbhd) | b > ncol(nbhd)){
          a <- 0
          b <- 0
        }
    }
    
    return(c(a,b))
  }
  
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

# returns the expected time to consensus over k square matrices
VM_func <- function(k){
  avg <- numeric(k)
  var <- numeric(k)
  sec <- numeric(k)
  
  for(i in 1:k){
    print(i)
    time <- VM(i, i, 100)
    avg[i] <- mean(time)
    var[i] <- var(time)
    sec[i] <- mean(time^2)
  }
  
  return(c(avg, var, sec))
}

# ----- #

# Confidence voter model - marginal
# two states are still 0 and 1
# there will be a matrix of same size that records each voters confidence to
# the corresponding location
CVM_marg <- function(row, col, s){
  
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
    
    return(c(a,b))
  }
  
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
CVM_extr <- function(row, col, s){
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
    
    return(c(a,b))
  }
  
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

# returns values for both marginal and extremal models
CVM_func <- function(fixed.row, col.range){
  avg_m <- numeric()
  avg_e <- numeric()
  var_m <- numeric()
  var_e <- numeric()
  sec_m <- numeric()
  sec_e <- numeric()
  
  for(i in col.range){
    print(i)
    time_m <- CVM_marg(fixed.row, i, 1000)
    time_e <- CVM_extr(fixed.row, i, 1000)
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

# Runs the same as VM, however this time we assume no boundaries
VM_noBound <- function(row,col,s){
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    
    index <- sample(1:4,1, replace = T)
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
    states <- sample(c(0,1), N, replace = TRUE)
    nbhd <- matrix(data = states, nrow = row, ncol = col)
    
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
      a <- lst[1]
      b <- lst[2]
      
      nbhd[i,j] <- nbhd[a,b]
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# input will be our fixed number of rows, and the various number of columns
# runs VM_noBound over the various sizes inputted by user
# returns avg, var, and second moment of each matrix size
VM_func_noBound <- function(fixed.row, col.range){
  avg <- numeric()
  var <- numeric()
  sec <- numeric()
  
  for(i in col.range){
    time <- VM_noBound(fixed.row, i, 1000)
    avg <- append(avg,mean(time))
    var <- append(var,var(time))
    sec <- append(sec,mean(time^2))
  }
  
  return(c(avg, var, sec))
}

# ----- #

# runs similiarly to CVM_marg, but
# state/confidence only change if neighbor is confident
CVM_marg_mod <- function(row,col,s){
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
    
    return(c(a,b))
  }
  
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
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# runs similiarly to CVM_extr, but
# state/confidence only change if neighbor is confident
CVM_extr_mod <- function(row, col, s){
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
    
    return(c(a,b))
  }
  
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
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# input will be our fixed number of rows, and the various number of columns
# runs CVM_marg_mod and CVM_extr_mod over the various sizes inputted by user
# returns avg, var, and second moment of each matrix size for each model
CVM_func_mod <- function(fixed.row, col.range){
  avg_m <- numeric()
  avg_e <- numeric()
  var_m <- numeric()
  var_e <- numeric()
  sec_m <- numeric()
  sec_e <- numeric()
  
  for(i in col.range){
    print(i)
    time_m <- CVM_marg_mod(fixed.row, i, 1000)
    time_e <- CVM_extr_mod(fixed.row, i, 1000)
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
VM_anim_special_1 <- function(row = 8, col = 12){
  
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
    
    return(c(a,b))
  }
  
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
VM_anim_special_2 <- function(row = 8, col = 12){
  
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
    
    return(c(a,b))
  }
  
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

# runs similarly to VM
# but returns a vector of which value wins consensus in each simulation
VM_win <- function(row,col,s, pOne = 0.5){
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    
    index <- sample(1:4,1, replace = T)
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
      lst <- neighbor(nbhd, i, j)
      a <- lst[1]
      b <- lst[2]
      
      nbhd[i,j] <- nbhd[a,b]
    
    }
  }
  print(paste("One wins at probability ", pOne,": ", sum(winner == 1)/s, sep = ""))
  return(winner)
}

# runs similarly to VM
# returns the percent of individuals exist with state '1'
VM_percents <- function(row,col, pOne = 0.5){
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    
    index <- sample(1:4,1, replace = T)
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
    lst <- neighbor(nbhd, i, j)
    a <- lst[1]
    b <- lst[2]
    
    nbhd[i,j] <- nbhd[a,b]
    
    percentOne <- append(percentOne, sum(nbhd == 1)/N)
  }
  
  return(percentOne)
}

# ----- #

# Individuals may now interact with neighbors "diagonal" to them
# Assumes no boundaries
VM_diag <- function(row, col, s, pOne = 0.5){
  
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0),
                      c(1,1),c(1,-1),c(-1,1),c(-1,-1))
    
    
    index <- sample(1:4,1, replace = T)
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

# input will be our fixed number of rows, and the various number of columns
# runs VM_diag over the various sizes inputted by user
# returns avg, var, and second moment of each matrix size
VM_diag_func <- function(fixed.row, col.range){
  avg <- numeric()
  var <- numeric()
  sec <- numeric()
  
  for(i in col.range){
    print(i)
    time <- VM_diag(fixed.row, i, 1000)
    avg <- append(avg,mean(time))
    var <- append(var,var(time))
    sec <- append(sec,mean(time^2))
  }
  
  return(c(avg, var, sec))
}

# row: number of rows
# col: number of columns
# s: number of observations
# pOne: probability that a spot in the matrix is initialized with '1'
# pChange.C: probability that an individual is willing to switch opinions after
         # interacting with a confident neighbor
# pChange.U: probability that an individual is willing to switch opinions after
         # interacting with an unsure neighbor
# We will assume boundaries, and that we cannot interact diagonally
my_conf_marg <- function(row, col, s, pOne = 0.5, pChange.C = 1, pChange.U = 0){
  # selects neighbor depending on if an individual is in the middle, edge, or corner
  neighbor <- function(nbhd,i,j){
    a <- 0
    b <- 0
    
    # possible movements (up/down/left/right)
    movements <- list(c(0,1),c(0,-1),c(-1,0),c(1,0))
    
    while(a <= 0 | b <= 0){
      index <- sample(1:4,1, replace = T)
      a <- i + movements[[index]][1]
      b <- j + movements[[index]][2]
      
      if( a > nrow(nbhd) | b > ncol(nbhd)){
        a <- 0
        b <- 0
      }
    }
    
    return(c(a,b))
  }
  
  N <- row*col
  consensusT <- numeric(s)
  
  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE, prob = c(1 - pOne, pOne))
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
      
      
      consensusT[k] <- consensusT[k] + 1
    }
  }
  
  return(consensusT)
}

# same parameters as my_conf_marg
my_conf_extr <- function(row, col, s, pOne = 0.5, pChange.C = 1, pChange.U = 0){
  
}