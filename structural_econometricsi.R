library(plyr)
library(stats4)
library(bbmle)

# Loading data
streco <- read.csv("C:/Users/benoit/Desktop/ENS/M2 MIE/S2/Structural eco/streco.csv", header=FALSE, stringsAsFactors=FALSE)
exp1 <- streco[-1,][1:6]
exp2 <- streco[-1,][7:12]
exp3 <- streco[-1,][13:18]
exp4 <- streco[-1,][19:24]
exp5 <- streco[-1,][25:30]
exp6 <- streco[-1,][31:36]
exp7 <- streco[-1,][37:42]
exp8 <- streco[-1,][43:48]

# Returns first time an individual chooses choice number 1 within an experiment
firstone <- function(v) {
  i = 1;
  while (v[i]< 1 & i < 7) {
    i = i+1;
  }
  return(i)
}

exp1bis <- apply(exp1, 1, function(x) firstone(x))
exp2bis <- apply(exp2, 1, function(x) firstone(x))
exp3bis <- apply(exp3, 1, function(x) firstone(x))
exp4bis <- apply(exp4, 1, function(x) firstone(x))
exp5bis <- apply(exp5, 1, function(x) firstone(x))
exp6bis <- apply(exp6, 1, function(x) firstone(x))
exp7bis <- apply(exp7, 1, function(x) firstone(x))
exp8bis <- apply(exp8, 1, function(x) firstone(x))

# Constructing the dataframe 
df <- data.frame(exp1bis, exp2bis, exp3bis, exp4bis, exp5bis, exp6bis, exp7bis)   

# Function to predict from a value of r on the 8th experiment
shouldChoose <- function(r) {
  i=1;
  plus=c(-50000000,5/100,1/10,2/10,5/10,1,2,5000000);
  while((1+plus[i]<(52+64*r)/(52+12*r))&i<7) {
    i=i+1;
  }
  return(i)
}

# Making predictions
results <- data.frame()
for (i in 2:381){ 
  tryCatch({
    likelihood <- function(beta, r, a) {
      x <- c(0,0,0,0,0,0,0)
      plus <- c(-50000000,5/100,1/10,2/10,5/10,1,2,5000000)
      x[1] <- log(pnorm(- (1/beta) + (1 / (1 + 4*r/52)) * ((4/52)*plus[df[i-1,1] + 1] + 1),0,a) - pnorm(- (1/beta) + (1 / (1 + 4*r/52)) * ((4/52)*plus[df[i-1,1]] + 1),0,a))
      x[2] <- log(pnorm(-(1/(1+1*r/52))+(1/(1+5*r/52))*((4/52)*plus[df[i-1,2]+1]+1),0,a)-pnorm(-(1/(1+1*r/52))+(1/(1+5*r/52))*((4/52)*plus[df[i-1,2]]+1),0,a))
      x[3] <- log(pnorm(-(1/(1+4*r/52))+(1/(1+8*r/52))*((4/52)*plus[df[i-1,3]+1]+1),0,a)-pnorm(-(1/(1+4*r/52))+(1/(1+8*r/52))*((4/52)*plus[df[i-1,3]]+1),0,a))
      x[4] <- log(pnorm(-(1/(1+12*r/52))+(1/(1+16*r/52))*((4/52)*plus[df[i-1,4]+1]+1),0,a)-pnorm(-(1/(1+12*r/52))+(1/(1+16*r/52))*((4/52)*plus[df[i-1,4]]+1),0,a))
      x[5] <- log(pnorm(- (1/beta)+(1/(1+r))*(plus[df[i-1,5]+1]+1),0,a)-pnorm(- (1/beta)+(1/(1+r))*(plus[df[i-1,5]]+1),0,a))
      x[6] <- log(pnorm(-(1/(1+r/52))+(1/(1+53*r/52))*(plus[df[i-1,6]+1]+1),0,a)-pnorm(-(1/(1+r/52))+(1/(1+53*r/52))*(plus[df[i-1,6]]+1),0,a))
      x[7] <- log(pnorm(-(1/(1+4*r/52))+(1/(1+56*r/52))*(plus[df[i-1,7]+1]+1),0,a)-pnorm(-(1/(1+4*r/52))+(1/(1+56*r/52))*(plus[df[i-1,7]]+1),0,a))
      -sum(x)
    }
    results[i-1,1] <- shouldChoose(mle2(minuslogl = likelihood, start = list(beta = 1, r = 0.1, a = 1))@coef['r'])
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  results[i-1,2] <- exp8bis[i-1]
}

# Counting successful predictions
success <- results[,1] - results[,2]
success <- success[!is.na(success)]
sum(success == 0)

# Counting times where prediction is really false (more than two apart)
sum(abs(success) > 2)

# Getting the values of beta, r and a for all individuals
parameter_values <- data.frame()
for (i in 2:381){ 
  tryCatch({
    likelihood <- function(beta, r, a) {
      x <- c(0,0,0,0,0,0,0)
      plus <- c(-50000000,5/100,1/10,2/10,5/10,1,2,5000000)
      x[1] <- log(pnorm(- (1/beta) + (1 / (1 + 4*r/52)) * ((4/52)*plus[df[i-1,1] + 1] + 1),0,a) - pnorm(- (1/beta) + (1 / (1 + 4*r/52)) * ((4/52)*plus[df[i-1,1]] + 1),0,a))
      x[2] <- log(pnorm(-(1/(1+1*r/52))+(1/(1+5*r/52))*((4/52)*plus[df[i-1,2]+1]+1),0,a)-pnorm(-(1/(1+1*r/52))+(1/(1+5*r/52))*((4/52)*plus[df[i-1,2]]+1),0,a))
      x[3] <- log(pnorm(-(1/(1+4*r/52))+(1/(1+8*r/52))*((4/52)*plus[df[i-1,3]+1]+1),0,a)-pnorm(-(1/(1+4*r/52))+(1/(1+8*r/52))*((4/52)*plus[df[i-1,3]]+1),0,a))
      x[4] <- log(pnorm(-(1/(1+12*r/52))+(1/(1+16*r/52))*((4/52)*plus[df[i-1,4]+1]+1),0,a)-pnorm(-(1/(1+12*r/52))+(1/(1+16*r/52))*((4/52)*plus[df[i-1,4]]+1),0,a))
      x[5] <- log(pnorm(- (1/beta)+(1/(1+r))*(plus[df[i-1,5]+1]+1),0,a)-pnorm(- (1/beta)+(1/(1+r))*(plus[df[i-1,5]]+1),0,a))
      x[6] <- log(pnorm(-(1/(1+r/52))+(1/(1+53*r/52))*(plus[df[i-1,6]+1]+1),0,a)-pnorm(-(1/(1+r/52))+(1/(1+53*r/52))*(plus[df[i-1,6]]+1),0,a))
      x[7] <- log(pnorm(-(1/(1+4*r/52))+(1/(1+56*r/52))*(plus[df[i-1,7]+1]+1),0,a)-pnorm(-(1/(1+4*r/52))+(1/(1+56*r/52))*(plus[df[i-1,7]]+1),0,a))
      -sum(x)
    }
    temp = mle2(minuslogl = likelihood, start = list(beta = 1, r = 0.1, a = 1))
    parameter_values[i-1,1] <- temp@coef['beta']
    parameter_values[i-1,2] <- temp@coef['r']
    parameter_values[i-1,3] <- temp@coef['a']
    confinttable = confint(temp)
    betaisnotone[i-1,1] <- (confinttable[1,1] > 1 || confinttable[1,2] < 1)
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
colnames(parameter_values) <- c('Beta', 'r', 'a')
# Results
parameter_values