library(purrr)
library(geepack)
library(MASS)
library(lme4)
library(geeM)


####parameters####
simu <- 200
I = 100 #number of participants
# I = rdunif(1, 10000, 2000)
T0 = 45 #number of follow-ups

sigma.a <- 0.1 # 0.5 
alpha <- rnorm(I, 0, sigma.a^2) #subject-specific intercept
beta.t = 0.130
beta.M = c(0, 4.53, 5.03, 6.02) 
beta.F = c(0, -1.99, -2.01, -2.84, -2.50) #how to specify?
true_para <- c(beta.t, beta.M[-1], beta.F[-1])
a <- 0.25 #1;25  dispersion parameter


####Simulate data####
Reward <- matrix(nrow = T0*I, ncol = simu)
Mot <- matrix(nrow = T0*I, ncol = simu)
Feed <- matrix(nrow = T0*I, ncol = simu)
t0 <- matrix(nrow = T0*I, ncol = simu)
id <- matrix(nrow = T0*I, ncol = simu)
for (s in 1:simu) {
  #Initial
  for (i in 1:I) {
    alpha.i <- alpha[i]
    for (t in 1:T0) {
      ####independent variable####
      Mot.ti <- rdunif(1, 4, 1) #Motivational message of user i at time t
      Feed.ti <- rdunif(1, 5, 1) #Feedback message of user i at time t
      Mot[(i-1)*T0 + t, s] <- Mot.ti
      Feed[(i-1)*T0 + t, s] <- Feed.ti
      t0[(i-1)*T0 + t, s] <- t
      id[(i-1)*T0 + t, s] <- i
      ####dependent variable####(Possion)
      lamb.ti <- exp(alpha.i + beta.t*t + beta.M[Mot.ti] + beta.F[Feed.ti])
      #theta.ti <- rgamma(1, shape = a, rate = 1/a) #overdispersion
      theta.ti <- 1 #without overdispersion
      Reward[(i-1)*T0 + t, s] <- rpois(1, lamb.ti*theta.ti)
    }
  }
}

####dataset####
Data <- list(Mot=Mot, Feed=Feed, Reward=Reward, id=id, t=t0, true=true_para)
####plot####


####Analysis####
Estresults <- function(Data, transformation) {
  PointEst <- function(Mot, Feed, t, Reward, id) {
    switch(transformation,
           Poi=coef(geeglm(Reward ~ t + factor(Mot) + factor(Feed), family = poisson, id = id, corstr = "independence"))[-1],
           Quasipoi=coef(geem(Reward ~ t + factor(Mot) + factor(Feed), family = quasipoisson, id = id, corstr = "independence"))[-1],
           NB=coef(geem(Reward ~ t + factor(Mot) + factor(Feed), family = negative.binomial(theta=2), id = id, corstr = "independence"))[-1],
           log=coef(lmer(log(Reward+1) ~ t + factor(Mot) + factor(Feed) + (1|id)))$id[1,-1],
           SR=coef(lmer(sqrt(Reward) ~ t + factor(Mot) + factor(Feed) + (1|id)))$id[1,-1],
           bcox=coef(lmer(((Reward^0.5 - 1)/0.5) ~ t + factor(Mot) + factor(Feed) + (1|id)))$id[1,-1])
  }
  Est <- matrix(nrow = 8, ncol = simu)
  Bias <- matrix(nrow = 8, ncol = simu)
  Relative.bias <- matrix(nrow = 8, ncol = simu)
  
  for (s in (1:simu)) {
    Est[,s]=unlist(PointEst(Data$Mot[,s], Data$Feed[,s], Data$t[,s], Data$Reward[,s], Data$id[,s]))
    Bias[,s]=Est[,s] - Data$true
    Relative.bias[,s]=Bias[,s]/Data$true
  }
  
  return(
  data.frame(
  true=Data$true,
  Meanest=apply(Est, 1, mean),
  EmpiricalVar=apply(Est, 1, var),
  MeanBias=apply(Bias, 1, mean),
  MeanRMSE=apply(Bias, 1, function(x) sqrt(mean(x^2))),
  MeanRelbias=apply(Relative.bias, 1, mean)
  )
  )
}

GetAnalyses=function(Data) {
  list(
    Data=Data,
    NB=Estresults(Data, transformation="NB"),
    Poi=Estresults(Data, transformation="Poi"),
    Quasipoi=Estresults(Data, transformation="Poi"),
    SR=Estresults(Data, transformation="SR"),
    log=Estresults(Data, transformation="log"),
    bcox=Estresults(Data, transformation="bcox")
  )
}

####Comparison results
results <- GetAnalyses(Data)

