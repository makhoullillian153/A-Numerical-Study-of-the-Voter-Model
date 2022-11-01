source("Voter model simulations.R")

time1 <- VM(5,5,100)
# result: 
# mean: 788.8189
# variance: 422908.7
# second moment: 1045102

fixed.row <- 8
col.range <- 1:15
k <- 15

lst <- VM_func(fixed.row, col.range)

avg <- lst[1:k]
var <- lst[(k+1):(2*k)]
sec <- lst[(2*k+1):(3*k)]

x <- c(1:k)^2

# Plots 1-3
{
  # average
  plot1 <- ggplot(mapping = aes(1:15,avg)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ I(x^2) + x) +
    xlab('Number of Columns') +
    ylab("Average Time to Absorption") +
    labs(title = "Fixed Number of Rows: 8")
  
  plot1
  
  x <- 1:15
  
  lm1 <- lm(avg ~ I(x^2) + x)
  
  # variance
  plot2 <- ggplot(mapping = aes(1:15, var)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,4,raw=TRUE)) +
    xlab('Number of Columns') +
    ylab("Variance of time to census") +
    labs(title = "Fixed Number of Rows: 8")
  
  plot2
  
  lm2 <- lm(var ~ poly(1:15,4,raw=TRUE))
  # y = 890500000
  # x0 should be close to 0; numbers are far too large
  
  # second moment
  plot3 <- ggplot(mapping = aes(x,sec)) + geom_point() +
    geom_smooth(method = "lm", formula = y ~ poly(x,4,raw=TRUE)) +
    xlab('Number of Columns') +
    ylab("Second moment of time to census") +
    labs(title = "Fixed Number of Rows: 8")
  
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
  plot4 <- ggplot(mapping = aes(x = 1:10)) + 
    geom_point(mapping = aes(y =avg_m, color = "Marginal")) +
    geom_point(mapping = aes(y = avg_e, color = "Extremal")) +
    #geom_smooth(mapping = aes(y = avg_m), method = "lm", formula = y ~ x, size = 0.5) +
    #geom_smooth(mapping = aes(y= avg_e), method = "lm", formula = y ~ x, color = "red", size = 0.5) +
    geom_line(mapping = aes(1:10, exp(1:10))) +
    xlab('Number of Columns') + 
    ylab("Average Time to Absorption") +
    scale_color_manual(name = "Legend",breaks = c("Marginal","Extremal"), 
                       values = c("Marginal" = "black", "Extremal" = "green"))
  
  plot4 
  # plots are pretty close but at N = 100 we can see the marginal distribution 
  # starting to diverge away from the extremal distribution
  
  lm4.1 <- lm(log(avg_m) ~ logx) # 126.463 - 27.664x + 2.537x^2
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

k <- 15
fixed.row <- 8
col.range <- 1:k

lst3 <- VM_func_noBound(fixed.row, col.range)

avg_mod <- lst3[1:k]
var_mod <- lst3[(k+1):(2*k)]
sec_mod <- lst3[(2*k+1):(3*k)]

x <- col.range

# Plots 7-9

{
  # average
  plot7 <- ggplot(mapping = aes(1:15, avg_mod)) + 
    geom_point() +
    #geom_point(mapping = aes(1:15, avg), color = "red")+
    geom_smooth(method = "lm", formula = y ~ x) +
    xlab('Number of columns') +
    ylab("Average time to census")
  
  plot7
  
  lm7 <- lm(avg_mod ~ c(1:15))
  
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

fixed.row <- 8
k <- 15
col.range <- 1:k
lst4 <- CVM_func_V2(fixed.row,col.range)

avg_m_mod <- lst4[1:k]
var_m_mod <- lst4[(k+1):(2*k)]
sec_m_mod <- lst4[(2*k+1):(3*k)]

avg_e_mod <- lst4[(3*k+1):(4*k)]
var_e_mod <- lst4[(4*k+1):(5*k)]
sec_e_mod <- lst4[(5*k+1):(6*k)]

# Plots 10-13

# compare the averages of original models to modified
{
  plot10 <- ggplot(mapping = aes(x = 1:10)) + 
    geom_point(mapping = aes(y = avg_m, color = "Original")) +
    geom_point(mapping = aes(y = avg_m_mod, color = "Version 2")) +
    #geom_smooth(mapping = aes(y = avg_m), method = "lm", formula = y ~ poly(x, 4, raw = TRUE), size = 0.5) +
    #geom_smooth(mapping = aes(y = avg_m_mod), method = "lm", formula = y ~ poly(x, 4, raw = TRUE), color = "red", size = 0.5) +
    xlab('Number of columns') + 
    ylab("Average time to consensus") +
    scale_color_manual(name = "Legend",breaks = c("Original","Version 2"), 
                       values = c("Original" = "black", "Version 2" = "green")) +
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
# mean: 450.586 (slightly than original version)
# second moment: 330141.6
# variance: 127241.1

fixed.row <- 8
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
plot17 <- ggplot(mapping = aes(x = 1:15)) + geom_point(mapping = aes(y = avg_diag, color = "Diagonal")) +
  geom_point(mapping = aes(y = avg_mod, color = "Classic")) +
  xlab('Number of columns') +
  ylab("Average time to census") +
  scale_color_manual(name = "Legend",breaks = c("Classic","Diagonal"), 
                     values = c("Classic" = "black", "Diagonal" = "red")) +
  labs(title = "Classic movements v/s diagonal",
       subtitle = "Assumes NO boundaries",
       caption = "Each simulation ran with 100000 observations")

plot17
# i really see no difference here

lm17 <- lm(avg_diag ~ poly(x,2,raw=TRUE))

lm17$coefficients
lm7$coefficients

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

lst6 <- couplingVM(6,10)

percent_classic <- lst6[1:(length(lst6)/2)]
percent_diag <- lst6[(length(lst6)/2 + 1):length(lst6)]

w <- which(percent_classic == percent_diag)

plot21 <- ggplot(mapping = aes(x = 1:(length(lst6)/2))) +
  geom_line(mapping = aes(y = percent_classic, color = "Classic")) +
  geom_line(mapping = aes(y = percent_diag, color = "Diagonal")) +
  geom_point(mapping = aes(x=w,y=percent_classic[w]), color = "blue", size = 2) +
  xlab("Time") +
  ylab("Percent of Particles labeled 'One'") +
  labs(title = "Comparing the Classic and Diagonal Models") + 
  scale_color_manual(name = "Legend",breaks = c("Classic","Diagonal"), 
                     values = c("Classic" = "black", "Diagonal" = "red"))

plot21

time6 <- my_conf_marg(5,5,1000)
# results:
  # mean: 6138.126 (takes much longer than without incorporating a unique probability of being convinced)
  # variance: 162074187
  # second moment: 199588704

fixed.row <- 8
k <- 15
col.range <- 1:k

lst7 <- my_conf_func(fixed.row,col.range)

avg_m_v3 <- lst7[1:k]
var_m_v3 <- lst7[(k+1):(2*k)]
sec_m_v3 <- lst7[(2*k+1):(3*k)]

avg_e_v3 <- lst7[(3*k+1):(4*k)]
var_e_v3 <- lst7[(4*k+1):(5*k)]
sec_e_v3 <- lst7[(5*k+1):(6*k)]

plot22 <- ggplot(mapping = aes(x = 1:15)) + 
  geom_point(mapping = aes(y = avg_m_mod, color = "1")) +
  geom_point(mapping = aes(y = avg_m_v3, color = "Varies")) +
  xlab('Number of columns') + 
  ylab("Average time to consensus") +
  scale_color_manual(name = "P(Cnvncd by C Neighbor)",breaks = c("1","Varies"), 
                     values = c("1" = "black", "Varies" = "green")) +
  ggtitle("Marginal: Rows = 8")

plot22

plot23 <- ggplot(mapping = aes(x = col.range)) + 
  geom_point(mapping = aes(y = avg_e_mod, color = "1")) +
  geom_point(mapping = aes(y = avg_e_v3, color = "Varies")) +
  xlab('Number of columns') + 
  ylab("Average time to consensus") +
  scale_color_manual(name = "Legend",breaks = c("1","Varies"), 
                     values = c("1" = "black", "Varies" = "green")) +
  ggtitle("Extremal: Rows = 8")

plot23

# fixed probability (0.5) of being convinced
lst8 <- CVM_func_mod(fixed.row,col.range, change = 0.5)

avg_m8 <- lst8[1:k]
var_m8 <- lst8[(k+1):(2*k)]
sec_m8 <- lst8[(2*k+1):(3*k)]

avg_e8 <- lst8[(3*k+1):(4*k)]
var_e8 <- lst8[(4*k+1):(5*k)]
sec_e8 <- lst8[(5*k+1):(6*k)]


# compare 3 different versions
plot24 <- ggplot(mapping = aes(x = 1:15)) + 
  geom_point(mapping = aes(y = avg_m_mod, color = "1")) +
  #geom_errorbar(aes(ymin = avg_m_mod-sqrt(var_m_mod), ymax= avg_m_mod+sqrt(var_m_mod), color = "1")) +
  geom_point(mapping = aes(y = avg_m_v3, color = "Varies")) +
  #geom_errorbar(aes(ymin = avg_m_v3-sqrt(var_m_v3), ymax= avg_m_v3+sqrt(var_m_v3), color = "Varies")) +
  geom_point(mapping = aes(y = avg_m8, color = "0.5")) +
  #geom_errorbar(aes(ymin = avg_m8-sqrt(var_m8), ymax= avg_m8+sqrt(var_m8), color = "0.5")) +
  xlab('Number of columns') + 
  ylab("Average time to consensus") +
  scale_color_manual(name = "P(Cnvncd by C Neighbor)", breaks = c("1","Varies","0.5"), 
                     values = c("1" = "black", "Varies" = "green", "0.5" = "red"))

plot24

plot25 <- ggplot(mapping = aes(x = 1:15)) + 
  geom_point(mapping = aes(y = avg_e_mod, color = "1")) +
  #geom_errorbar(aes(ymin = avg_e_mod-sqrt(var_e_mod), ymax= avg_e_mod+sqrt(var_e_mod), color = "1")) +
  geom_point(mapping = aes(y = avg_e_v3, color = "Varies")) +
  #geom_errorbar(aes(ymin = avg_e_v3-sqrt(var_e_v3), ymax= avg_e_v3+sqrt(var_e_v3), color = "Varies")) +
  geom_point(mapping = aes(y = avg_e8, color = "0.5")) +
  #geom_errorbar(aes(ymin = avg_e8-sqrt(var_e8), ymax= avg_e8+sqrt(var_e8), color = "0.5")) +
  xlab('Number of columns') + 
  ylab("Average time to consensus") +
  scale_color_manual(name = "P(Cnvncd by C Neighbor)", breaks = c("1","Varies","0.5"), 
                     values = c("1" = "black", "Varies" = "green", "0.5" = "red"))

plot25

# now lets look at the coefficients of each model
# these values really aren't that reliable with such a low number of observations

# linear section - marginal, extremal, respectively
x <- 1:fixed.row

# 0 < P(Convinced) < 1
lm20.1 <- lm(avg_m_v3[1:fixed.row] ~ x) # -5152 + 3466x
lm20.2 <- lm(avg_e_v3[1:fixed.row] ~ x) # -8332 + 4540x

# P(Convinced) = 0.5
lm21.1 <- lm(avg_m8[1:fixed.row] ~ x) # -1304 + 1224x
lm21.2 <- lm(avg_e8[1:fixed.row] ~ x) # -928 + 763.7x

# P(Convinced) = 1
lm22.1 <- lm(avg_m_mod[1:fixed.row] ~ x) # -655.8 + 605.4x
lm22.2 <- lm(avg_e_mod[1:fixed.row] ~ x) # -452.2 + 376x

# comparing SLOPES
# Coeff21/coeff22 -> 2/1 (very close to 2, actually)
# interpretation:
  # the time to consensus in model 'mod' grows twice as fast as the time to consensus
  # from the original confidence model, for both marginal and extremal

# Coeff20.1/21.1 -> 3/1
# Coeff20.2/21.2 -> 6/1

# Coeff20.1/Coeff22.1 -> 6/1
# Coeff20.2/Coeff22.2 -> 12/1

# it looks like the slope ratios for extremal is double that of marginal for 
# model 20 to 21, and 20 to 22. i wonder how much closer/further these values will
# be from their estimates as i increase the number of observations.

# "exponential" section - marginal, extremal, respectively
x <- (fixed.row+1):k

# no compelling results from this analysis, the more compelling results stem from 
# the linear portion of the results

## ignoring intercept, not our main interest when comparing these coefficients

# 0 < P(Convinced) < 1
lm20.3 <- lm(avg_m_v3[x] ~ poly(x,2,raw = TRUE)) # -44263x + 3658x^2
lm20.4 <- lm(avg_e_v3[x] ~ poly(x,2,raw = TRUE)) # -131556x + 9181x^2

# P(Convinced) = 0.5
lm21.3 <- lm(avg_m8[x] ~ poly(x,2,raw = TRUE)) # -17074x + 1493x^2
lm21.4 <- lm(avg_e8[x] ~ poly(x,2,raw = TRUE)) # -2876.1x + 400.5x^2

# P(Convinced) = 1
lm22.3 <- lm(avg_m_mod[x] ~ poly(x,2,raw = TRUE)) # -12764x + 1007x^2
lm22.4 <- lm(avg_e_mod[x] ~ poly(x,2,raw = TRUE)) # -2050.5x + 222.1x^2

# comparing X1, X2
# Coeff21.3/coeff22.3 -> 4/3 , 3/2
# Coeff21.4/coeff22.4 -> 1.4 , 18/10

# Coeff20.3/21.3 -> 2.592, 2.45 (ratios converging to 5/2?)
# Coeff20.4/21.4 -> 45.741, 22.92385

# Coeff20.3/Coeff22.3 -> 3.4678, 2.530288
# Coeff20.4/Coeff22.4 -> 64.158, 41.337

time7 <- noConsensus_percent(3,2)

plot26 <- ggplot(mapping = aes(x = 1:length(time7), y = time7)) + 
  geom_line() + 
  geom_hline(mapping = aes(yintercept = mean(time7)),color = 'red') +
  xlab("Time") +
  ylab("Percent of Particles labeled 'One'") + 
  ylim(0,1)

set.seed(1)
{time9 <- noConsensus_time(3,3,1000)
plot27 <- ggplot(mapping = aes(x = 1:length(time9), y = time9)) + 
  geom_line() + 
  geom_hline(mapping = aes(yintercept = mean(time9)), color = "red") +
  xlab("Time") +
  ylab("Percent of Particles labeled 'One'") +
  ylim(0,1)
plot27

}

plot28 <- ggplot(mapping = aes(x = c(51,106.1,186.56,5.9,18.1,42.4,133.4,302.4223),
                               y = c(1300.5,10482.2,95457.29,26.1256,160,795,16897.84,))) +
  geom_point()

mean(VM(3,3,10000)^2) # ~51
mean(V4_time(3,3,1000)) # 1263.61
#FOR 3x3
# 51*51*(1/2) = 1300.5

mean(VM(3,4,20000)) # 106.1
mean(V4_time(3,4,500)) # 9666.416
# 106*106*(2/2) = 11,236

mean(VM(3,5,5000)) # 186.56
mean(V4_time(3,5,100)) # 95457.29
# 186.6*186.6*(6/2) = 104,458

mean(VM(3,6,10000)) # 302.4223
mean(V4_time(3,6,500)) # 722845.9
# 302.42*302.42*(24/2) = 1,097,494.28

95457.29/(302.4*302.4)

mean(VM(2,2,10000)) # 5.9
mean(V4_time(2,2,10000)) # 26.1256
# 5.9*5.9*(3/4) = 26.1075

mean(VM(3,2,10000)) # 18.1
mean(V4_time(3,2,10000)) # 160
#18.1*18.1*(1/2) = 163.805

mean(VM(2,1,10000)) # 0.5
mean(V4_time(2,1,10000)) # 1
# 0.5*0.5*4 = 1

mean(VM(2,4,10000)) # 42.4
mean(V4_time(2,4,5000)) # 795

mean(VM(2,6,10000)) # 133.4
mean(V4_time(2,6,50)) # 16897.84
