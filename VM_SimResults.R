source("Voter model simulations.R")

time1 <- VM(5,5,100)
# result: 
# mean: 788.8189
# variance: 422908.7
# second moment: 1045102

k <- 10

lst <- VM_func(k)

avg <- lst[1:k]
var <- lst[(k+1):(2*k)]
sec <- lst[(2*k+1):(3*k)]

x <- c(1:k)^2

# Plots 1-3
{
  # average
  plot1 <- ggplot(mapping = aes(x,avg)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,2,raw=TRUE)) +
    xlab('Number of "residents"') +
    ylab("Average time to census")
  
  plot1
  
  lm1 <- lm(avg ~ poly(x,2,raw=TRUE))
  # y = 102.004 - 22.624x +1.893x^2
  
  # variance
  plot2 <- ggplot(mapping = aes(x, var)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,4,raw=TRUE)) +
    xlab('Number of "Residents"') +
    ylab("Variance of time to census")
  
  plot2
  
  lm2 <- lm(var ~ poly(x,4,raw=TRUE))
  # y = 890500000
  # x0 should be close to 0; numbers are far too large
  
  # second moment
  plot3 <- ggplot(mapping = aes(x,sec)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,4,raw=TRUE)) +
    xlab('Number of "Residents"') +
    ylab("Second moment of time to census")
  
  plot3
  
  lm3 <- lm(sec ~ poly(x, 4, raw = TRUE))
  # same issue: numbers far too large
  
  # log-log plots are all relatively linear
}

time2 <- CVM_marg(5,5,1000)
# result:
# mean: 1142.872
# variance: 1085055
# second moment: 2391102

time3 <- CVM_extr(5,5,1000)
# results:
# mean: 1185.004
# variance: 1164857
# second moment: 2568976
# values are slightly larger; surprising because unsure individuals jump straight
# to confident if they meet a differing opinion

k <- 10
fixed.row <- 6
col.range <- 1:k

lst2 <- CVM_func(fixed.row, col.range)

avg_m <- lst2[1:k]
var_m <- lst2[(k+1):(2*k)]
sec_m <- lst2[(2*k+1):(3*k)]

avg_e <- lst2[(3*k+1):(4*k)]
var_e <- lst2[(4*k+1):(5*k)]
sec_e <- lst2[(5*k+1):(6*k)]


# Plots 4-6

# let black dots represent marginal and green dots represent extremal
{
  x <- c(1:k)^2
  
  # average
  plot4 <- ggplot() + geom_point(mapping = aes(x, avg_m, color = "Marginal")) +
    geom_point(mapping = aes(x, avg_e, color = "Extremal")) +
    geom_smooth(mapping = aes(x, avg_m), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), size = 0.5) +
    geom_smooth(mapping = aes(x, avg_e), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "red", size = 0.5) +
    xlab('Number of "residents"') + 
    ylab("Average time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Marginal","Extremal"), 
                       values = c("Marginal" = "black", "Extremal" = "green"))
  
  plot4 
  # plots are pretty close but at N = 100 we can see the marginal distribution 
  # starting to diverge away from the extremal distribution
  
  lm4.1 <- lm(avg_m ~ poly(x,2, raw = TRUE)) # 126.463 - 27.664x + 2.537x^2
  lm4.2 <- lm(avg_e ~ poly(x,2, raw = TRUE)) # 84.502 - 21.838x + 2.448x^2
  
  # variance
  plot5 <- ggplot() + geom_point(mapping = aes(x, var_m, color = "Marginal")) +
    geom_point(mapping = aes(x, var_e, color = "Extremal")) +
    geom_smooth(mapping = aes(x, var_m), method = "lm", formula = y ~ poly(x, 4, raw = TRUE)) +
    geom_smooth(mapping = aes(x, var_e), method = "lm", formula = y ~ poly(x, 4, raw = TRUE), color = "red") +
    xlab('Number of "residents"') + 
    ylab("Variance of time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Marginal","Extremal"), 
                       values = c("Marginal" = "black", "Extremal" = "green"))
  
  plot5
  
  # numbers remain too large for variance
  lm5.1 <- lm(var_m ~ poly(x,4, raw = TRUE)) # 
  lm5.2 <- lm(var_e ~ poly(x,4, raw = TRUE)) # 
  
  # second moment
  plot6 <- ggplot() + geom_point(mapping = aes(x, sec_m, color = "Marginal")) +
    geom_point(mapping = aes(x, sec_e, color = "Extremal")) +
    geom_smooth(mapping = aes(x, sec_m), method = "lm", formula = y ~ poly(x, 4, raw = TRUE)) +
    geom_smooth(mapping = aes(x, sec_e), method = "lm", formula = y ~ poly(x, 4, raw = TRUE), color = "red") +
    xlab('Number of "residents"') + 
    ylab("Second moment of time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Marginal","Extremal"), 
                       values = c("Marginal" = "black", "Extremal" = "green"))
  
  plot6
  
  # numbers remain too large for second moment
  lm6.1 <- lm(sec_m ~ poly(x,4, raw = TRUE)) # 
  lm6.2 <- lm(sec_e ~ poly(x,4, raw = TRUE)) # 
}

time4 <- VM_noBound(2,2,1000)
# results:
# mean: 588.1 - faster
# variance: 160875.8
# second moment: 505128.6

fixed.row <- 6
col.range <- 1:10

lst3 <- VM_func_noBound(fixed.row, col.range)

avg_mod <- lst3[1:10]
var_mod <- lst3[(10+1):(2*10)]
sec_mod <- lst3[(2*10+1):(3*10)]

x <- c(1:10)

# Plots 7-9

{
  # average
  plot7 <- ggplot(mapping = aes(x,avg_mod)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,2,raw=TRUE)) +
    xlab('Number of columns') +
    ylab("Average time to census")
  
  plot7
  
  lm7 <- lm(avg_mod ~ poly(x,2,raw=TRUE))
  
  # variance
  plot8 <- ggplot(mapping = aes(x, var_mod)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,3,raw=TRUE)) +
    xlab('Number of columns') +
    ylab("Variance of time to census")
  
  plot8
  # variance at col = 9 is  higher than the others
  # col = 8 is lower
  
  lm8 <- lm(var_mod ~ poly(x,3,raw=TRUE))
  
  # second moment
  plot9 <- ggplot(mapping = aes(x,sec_mod)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,3,raw=TRUE)) +
    xlab('Number of columns') +
    ylab("Second moment of time to census")
  
  plot9
  
  lm9 <- lm(sec_mod ~ poly(x, 3, raw = TRUE))
}

time5 <- CVM_marg_mod(5,5,1000)
# results
# mean: 1414.98
# takes slightly longer than original marginal model (about 200 more steps, 20%)

fixed.row <- 6
k <- 10
col.range <- 1:k
lst4 <- CVM_func_mod(fixed.row,col.range)

avg_m_mod <- lst4[1:k]
var_m_mod <- lst4[(k+1):(2*k)]
sec_m_mod <- lst4[(2*k+1):(3*k)]

avg_e_mod <- lst4[(3*k+1):(4*k)]
var_e_mod <- lst4[(4*k+1):(5*k)]
sec_e_mod <- lst4[(5*k+1):(6*k)]

# Plots 10-13

# compare the averages of original models to modified
{
  plot10 <- ggplot(mapping = aes(x = col.range)) + 
    geom_point(mapping = aes(y = avg_m, color = "Original")) +
    geom_point(mapping = aes(y = avg_m_mod, color = "Modified")) +
    #geom_smooth(mapping = aes(y = avg_m), method = "lm", formula = y ~ poly(x, 4, raw = TRUE), size = 0.5) +
    #geom_smooth(mapping = aes(y = avg_m_mod), method = "lm", formula = y ~ poly(x, 4, raw = TRUE), color = "red", size = 0.5) +
    xlab('Number of columns') + 
    ylab("Average time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Original","Modified"), 
                       values = c("Original" = "black", "Modified" = "green")) +
    ggtitle("Marginal: Rows = 6")
  
  plot10
  
  plot11 <- ggplot(mapping = aes(x = col.range)) + 
    geom_point(mapping = aes(y = avg_e, color = "Original")) +
    geom_point(mapping = aes(y = avg_e_mod, color = "Modified")) +
    #geom_smooth(mapping = aes(y = avg_e), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), size = 0.5) +
    #geom_smooth(mapping = aes(y = avg_e_mod), method = "lm", formula = y ~ poly(x, 2, raw = TRUE), color = "red", size = 0.5) +
    xlab('Number of columns') + 
    ylab("Average time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Original","Modified"), 
                       values = c("Original" = "black", "Modified" = "green")) +
    ggtitle("Extremal: Rows = 6")
  
  plot11 # trends are not as consistent with the marginal model. maybe my 
  # conditioning in the simulations were wrong?
  
  # may be because voters jump straight to confident so we easily have
  # a growing confident population, making the average times similar
  
  plot12 <- ggplot(mapping = aes(x = col.range)) + 
    geom_point(mapping = aes(y = avg_e_mod, color = "Extremal")) +
    geom_point(mapping = aes(y = avg_m_mod, color = "Marginal")) +
    xlab('Number of columns') + 
    ylab("Average time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Extremal","Marginal"), 
                       values = c("Extremal" = "black", "Marginal" = "green")) +
    ggtitle("Modified: Rows = 6")
  
  plot12 # trends consistent with original: marginal values > extremal
  
  plot13 <- ggplot(mapping = aes(x = col.range)) + 
    geom_point(mapping = aes(y = avg_e, color = "Extremal")) +
    geom_point(mapping = aes(y = avg_m, color = "Marginal")) +
    xlab('Number of columns') + 
    ylab("Average time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Extremal","Marginal"), 
                       values = c("Extremal" = "black", "Marginal" = "green")) +
    ggtitle("Original: Rows = 6")
  
  plot13
}

# Plots 14-16 (VM_percents)

{
p <- 0.5
percents <- VM_percents(8,12, p)

plot14 <- ggplot() + 
  geom_line(mapping = aes(x = 1:length(percents), y = percents)) + 
  geom_hline(yintercept = p, color = 'red') +
  xlab("index") +
  ylab("Percent of Particles labeled 'One'")

plot14

p <- 0.75

percents <- VM_percents(8,12, p)

plot15 <- ggplot() + 
  geom_line(mapping = aes(x = 1:length(percents), y = percents)) + 
  geom_hline(yintercept = p, color = 'red') +
  xlab("index") +
  ylab("Percent of Particles labeled 'One'")

plot15

p <- 0.25

percents <- VM_percents(8,12, p)

plot16 <- ggplot() + 
  geom_line(mapping = aes(x = 1:length(percents), y = percents)) + 
  geom_hline(yintercept = p, color = 'red') +
  xlab("index") +
  ylab("Percent of Particles labeled 'One'")

plot16
}

time_diag <- VM_diag(5,5,1000)
# results
# mean: 497.211 (much faster than original version)
# second moment: 132841
# variance: 337938

fixed.row <- 6
k <- 15
col.range <- 1:k
lst5 <- VM_diag_func(fixed.row, col.range)
avg_diag <- lst5[1:k]
var_diag <- lst5[(k+1):(2*k)]
sec_diag <- lst5[(2*k+1):(3*k)]

x <- col.range

# Plots 17-19

{
# average
plot17 <- ggplot(mapping = aes(x = x)) + geom_point(mapping = aes(y = avg_diag, color = "Diagonal")) +
  geom_point(mapping = aes(y = avg_mod, color = "Classic")) +
  xlab('Number of columns') +
  ylab("Average time to census") +
  scale_color_manual(name = "Legend",breaks = c("Classic","Diagonal"), 
                     values = c("Classic" = "black", "Diagonal" = "red")) +
  labs(title = "Classic movements, v/s diagonal",
       subtitle = "Assumes NO boundaries",
       caption = "Each simulation ran with 1000 observations")

plot17

lm17 <- lm(avg_diag ~ poly(x,2,raw=TRUE))

# variance
plot18 <- ggplot(mapping = aes(x, var_diag)) + geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x,3,raw=TRUE)) +
  xlab('Number of columns') +
  ylab("Variance of time to census") +
  ggtitle("Diagonal interaction allowed")

plot18
# variance at col = 9 is  higher than the others
# col = 8 is lower

lm18 <- lm(var_diag ~ poly(x,3,raw=TRUE))

# second moment
plot19 <- ggplot(mapping = aes(x,sec_diag)) + geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x,3,raw=TRUE)) +
  xlab('Number of columns') +
  ylab("Second moment of time to census") +
  ggtitle("Diagonal interaction allowed")

plot19

lm19 <- lm(sec_diag ~ poly(x, 3, raw = TRUE))

# residual plot of plot17
plot20 <- ggplot(mapping = aes(x = x, y = avg_diag - avg_mod)) +
  geom_point()
}