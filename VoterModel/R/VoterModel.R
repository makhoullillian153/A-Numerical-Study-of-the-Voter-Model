#' @title neighbor
#' @description Uniformly selects neighbor of an individual using classic movements
#' @param nbhd Matrix representing the different individuals
#' @param i Row number of individual
#' @param j Column number of individual
#' @param boundaries If TRUE, we assume boundaries, otherwise edges and corners are linked to the other side of the matrix. Default is set to TRUE.
#' @return Row and column number of a randomly selected neighbor
#' @export
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
  }


  return(c(a,b))
}

#' @title Classic voter model
#' @description Used in VM theory.R file. Assumes boundaries.
#' @param row Number of rows
#' @param col Number of columns
#' @param s Number of observations
#' @param pOne Probability that a spot in the matrix is initialized with '1'. Default is 0.5.
#' @return Time to consensus of each observation
#' @export
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

      if(nbhd[i,j] != nbhd[a,b]){
        nbhd[i,j] <- nbhd[a,b]
        consensusT[k] <- consensusT[k] + 1
      }
    }
  }

  return(consensusT)
}

#' @title VM_func
#' @description Models the time to consensus of VM as a function over a number of columns
#' @param fixed.row Number of rows that will remain fixed through each simulation
#' @param col.range Range of columns
#' @return Average of, variance of, second moment of consensus of each dimension (3 vectors total)
#' @export
VM_func <- function(fixed.row, col.range){
  avg <- numeric()
  var <- numeric()
  sec <- numeric()
  j <- 1
  for(i in col.range){
    print(i)
    time <- VM(fixed.row, i, 10000)
    avg[j] <- mean(time)
    var[j] <- var(time)
    sec[j] <- mean(time^2)
    j <- j + 1
  }

  return(c(avg, var, sec))
}

# ----- #

#' @title Confidence voter model - marginal
#' @description Assumes boundaries
#' @param row number of rows
#' @param col number of columns
#' @param s number of observations
#' @return Time to consensus of each observation
#' @export
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

#' @title Confidence voter model - extremal
#' @description Assumes boundaries
#' @param row number of rows
#' @param col number of columns
#' @param s number of observations
#' @return Time to consensus of each observation
#' @export
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

#' @title CVM_func
#' @description Models the time to consensus of the confidence voter models as a function over a number of columns
#' @param fixed.row Number of rows that will remain fixed through each simulation
#' @param col.range Range of columns. Each will be used twice: once in marginal, and once in extremal case.
#' @return Average of, variance of, second moment of consensus of each dimension, for both marginal and extremal cases (6 vectors total)
#' @export
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

#' @title Confidence voter model - marginal (Version 2)
#' @description Runs similarly to CVM_marg, but state/confidence only change if neighbor is confident. We can also adjust the probability that an individual is convinced by a neighbor (conditioned on that neighbor being confident). Assumes boundaries
#' @param row Number of rows
#' @param col Number of columns
#' @param s Number of observations
#' @param change Probability that an individual will be convinced by a confident neighbor (constant for each individual). Default is 1
#' @return Time to consensus of each observation
#' @export
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

#' @title Confidence voter model - extremal (Version 2)
#' @description Runs similarly to CVM_extr, but state/confidence only change if neighbor is confident. We can also adjust the probability that an individual is convinced by a neighbor (conditioned on that neighbor being confident). Assumes boundaries
#' @param row Number of rows
#' @param col Number of columns
#' @param s Number of observations
#' @param change Probability that an individual will be convinced by a confident neighbor (constant for each individual). Default is 1
#' @return Time to consensus of each observation
#' @export
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

#' @title CVM_func_V2
#' @description Models the time to consensus of the confidence voter models (Version 2) as a function over a number of column
#' @param fixed.row Number of rows that will remain fixed through each simulation
#' @param col.range Range of columns. Each will be used twice: once in marginal, and once in extremal case.
#' @param change Probability that an individual will be convinced by a confident neighbor (constant for each individual). Default is 1
#' @param returns Average of, variance of, second moment of consensus of each dimension, for both marginal and extremal cases (6 vectors total)
#' @export
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

#' @title Confidence voter model - marginal (Version 3)
#' @description Each individual will have their own probability of being influenced by a confident voter. Assumes boundaries, and only classic movements
#' @param row Number of rows
#' @param col Number of columns
#' @param s Number of observations
#' @param pOne Probability that a spot in the matrix is initialized with '1'. Default is 0.5
#' @return Time to consensus of each observation
#' @export
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

#' @title Confidence voter model - extremal (Version 3)
#' @description Each individual will have their own probability of being influenced by a confident voter. Assumes boundaries, and only classic movements
#' @param row Number of rows
#' @param col Number of columns
#' @param s Number of observations
#' @param pOne Probability that a spot in the matrix is initialized with '1'. Default is 0.5
#' @return Time to consensus of each observation
#' @export
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

#' @title CVM_func_V3
#' @description Models the time to consensus of the confidence voter models (Version 3) as a function over a number of column
#' @param fixed.row Number of rows that will remain fixed through each simulation
#' @param col.range Range of columns. Each will be used twice: once in marginal, and once in extremal case.
#' @return: Average of, variance of, second moment of consensus of each dimension, for both marginal and extremal cases (6 vectors total)
#' @export
CVM_func_V3 <- function(fixed.row, col.range){
  avg_m <- numeric()
  avg_e <- numeric()
  var_m <- numeric()
  var_e <- numeric()
  sec_m <- numeric()
  sec_e <- numeric()

  for(i in col.range){
    print(i)
    time_m <- CVM_marg_V3(fixed.row, i, 2500)
    time_e <- CVM_extr_V3(fixed.row, i, 2500)
    avg_m[i] <- mean(time_m)
    avg_e[i] <- mean(time_e)

    var_m[i] <- var(time_m)
    var_e[i] <- var(time_e)

    sec_m[i] <- mean(time_m^2)
    sec_e[i] <- mean(time_e^2)
  }



  return(c(avg_m, var_m, sec_m, avg_e, var_e, sec_e))
}

#' @title My modified voter model, percents
#' @description Used in VM theory.R file. Here, the probability of being convinced is dependent on the proportion of agreeing voters. That is, if we have 3 individuals labeled "0" and 7 voters labeled "1", a "0" voter has a probability of 3/10 of being convinced, and a "1" voter has a probability of 7/10 of being convinced. We only run 1 observation. Assumes boundaries
#' @param row Number of rows
#' @param col Number of columns
#' @param pOne Probability that a spot in the matrix is initialized with '1'. Default is 0.5
#' @return: Percent of population labeled "1" throughout the observation
#' @export
myModified_percent <- function(row, col){

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


#' @title My modified voter model, times
#' @description used in VM theory.R file. Here, the probability of being convinced is dependent on the proportion of agreeing voters. That is, if we have 3 individuals labeled "0" and 7 voters labeled "1", a "0" voter has a probability of 3/10 of being convinced, and a "1" voter has a probability of 7/10 of being convinced. Assumes boundaries
#' @param row Number of rows
#' @param col Number of columns
#' @param s Number of observations
#' @param pOne Probability that a spot in the matrix is initialized with '1'. Default is 0.5
#' @return Time to consensus of each observation
#' @export
myModified_time <- function(row, col, s, pOne = 0.5){
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

#' @title Classic model on a complete graph
#' @description Looks at the classic voter model, but on a complete graph, rather than a square lattice graph.
#' @param N Number of individuals in population
#' @param s Number of observations
#' @param pOne Probability that a spot in the matrix is initialized with '1'. Default is 0.5
#' @return Time to consensus of each observation
#' @export
VM_complete <- function(N, s, pOne = 0.5){
  consensusT <- numeric(s)

  for(k in 1:s){
    states <- sample(c(0,1), N, replace = TRUE, prob = c(1-pOne, pOne))
    nbhd <- matrix(data = states, nrow = 1, ncol = N)

    while(TRUE){

      if(sum(nbhd) == N | sum(nbhd) == 0){
        break
      }

      # randomly select an individual
      i <- sample(1:N, 1)

      # select a neighbor - complete graph
      n <- sample(c(1:N)[-i], 1)

      nbhd[1,i] <- nbhd[1,n]
      consensusT[k] <- consensusT[k] + 1
    }
  }

  return(consensusT)
}

#' @title VM_complete_func
#' @description Models the time to consensus of VM_complete as a function over the initial density of agreeing voters.
#' @param row Number of rows
#' @param col Number of columns
#' @param s Number of observations
#' @return Average time to consensus of each initial density
#' @export
VM_complete_func <- function(row, col, s){
  emp <- numeric()
  N <- row*col

  for(i in 0:N){
    print(i)
    emp <- append(emp, mean(VM_complete(row, col, s, pOne = i/N)))
  }
  return(emp)
}
