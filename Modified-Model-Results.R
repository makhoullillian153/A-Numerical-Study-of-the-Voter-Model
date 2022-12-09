source("Voter-model-simulations.R")

## Minimum dimension 2

mean(VM(2,1,10000)) # 0.5
mean(myModified_time(2,1,10000)) # 1
# 0.5*0.5*4 = 1

# mean(VM_complete(2,2,10000)) # 3.02
# mean(myModified_complete(2,2,10000)) # 17.6967

mean(VM(2,2,10000)) # 5.9
mean(myModified_time(2,2,1000)) # 26.1256
# 5.9*5.9*(3/4) = 26.1075

# mean(VM_complete(2,3,10000)) # 7.4383
# mean(myModified_complete(2,3,1000)) # 102.37
mean(VM(2,3,10000)) # 18.1
mean(myModified_time(3,2,10000)) # 160
#18.1*18.1*(1/2) = 163.805

(18.1*18.1/2)-18.1

mean(VM(2,4,10000)) # 42.4
mean(myModified_time(2,4,5000)) # 795
# 42.4*42.4/2 = 898.88

42.4*42.4/2 - 42.4

mean(VM(2,5,10000)) # 78.4233
mean(myModified_time(2,5,5000)) # 3681.061
# 78.4*78.4/2 = 3073.28

78.4*78.4/2 - 78.4

mean(VM(2,6,10000)) # 133.4
mean(myModified_time(2,6,300)) # 16897.84
# 133.4*133.4 = 17795.6

(133.4*133.4) - 133.4

## Minimum dimension 3

# mean(VM_complete(3,3,10000)) # 47.3622
# mean(myModified_complete(3,3,1000)) # 1370.486

mean(VM(3,3,10000)) # 50.9
mean(myModified_time(3,3,1000)) # 1263.61
# 51*51*(1/2) - 51= 1249.5

51*51/2 - 51 # 1249.5

mean(VM(3,4,20000)) # 106.1
mean(myModified_time(3,4,500)) # 9666.416
# 106*106*(2/2) - 106*12 = 11,130

106*106*(2/2) - 106*(4^2)

mean(VM(3,5,10000)) # 186.56
mean(myModified_time(3,5,300)) # 83756.84
# 186.6*186.6*(6/2) = 104,458

186.6*186.6*(6/2) - 186.56*(5^3)

mean(VM(3,6,10000)) # 302.4223
mean(myModified_time(3,6,500)) # 722845.9
# 302.42*302.42*(24/2) = 1,097,494.28

302.42*302.42*(24/2)-(302)*(6^(4))

# mean(VM_complete(3,7,10000))
mean(VM(3,7,10000)) # 454.0221
mean(myModified_time(3,7,750)) # 6100133

## Minimum dimension 4

mean(VM(4,4,10000)) # 206.976
mean(myModified_time(4,4,1000)) # 144260

mean(VM(4,5,10000)) # 364.2978
mean(myModified_time(4,5,1000)) # 2272649

mean(VM(4,6,10000)) # 560.2132

mean(VM(4,7,10000)) # 826.2581

## Minimum dimension 5

mean(VM(5,5,10000)) # 606.671

mean(VM(5,6,10000)) # 937.198

df <- read.csv("Model data - Sheet1.csv")
VM.classic <- df$Classic.Time
VM.modified <- df$Modified.Time
VM.classic[11] <- 364.2978
VM.modified[9] <- 6100133

plot(VM.classic,VM.modified)
plot(VM.classic,log(VM.modified))
points(VM.classic, 6.154987 + 0.02335*VM.classic, col = "red")

plot(log(VM.classic),log(log(VM.modified)))

plot32 <- ggplot(mapping = aes(VM.classic, VM.modified)) +
  geom_point() +
  xlab(expression('T'[C])) +
  ylab(expression('T'[M]))

plot32

lm(log(VM.modified[1:4]) ~ log(VM.classic[1:4])) #slope: 2.322
lm(log(VM.modified[5:9]) ~ log(VM.classic[5:9])) #slope: 3.870
lm(log(VM.modified[10:11]) ~ log(VM.classic[10:11])) #slope: 4.775

plot29 <- ggplot(mapping = aes(x = log(VM.classic[1:4]), y = log(VM.modified[1:4]))) +
  geom_point() + 
  xlab(expression('log(T'[C]*')')) +
  ylab(expression('log(T'[M]*')'))

plot30 <- ggplot(mapping = aes(x = log(VM.classic[5:9]), y = log(VM.modified[5:9]))) +
  geom_point() + 
  xlab(expression('log(T'[C]*')')) +
  ylab(expression('log(T'[M]*')'))

plot31 <- ggplot(mapping = aes(x = log(VM.classic[10:11]), y = log(VM.modified[10:11]))) +
  geom_point() + 
  xlab(expression('log(T'[C]*')')) +
  ylab(expression('log(T'[M]*')'))
