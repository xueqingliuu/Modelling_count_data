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

sigma.a <- 0.5 #0.5
alpha <- rnorm(I, 5.4, sigma.a^2) #subject-specific intercept
beta.t = -0.0268 
beta.M = c(0, 2.81, 5.70, 4.21) 
beta.F = c(0, 1.07, 1.03, -3.49, -3.01) #how to specify?
true_para <- c(beta.t, beta.M[-1], beta.F[-1])
#a <- 10 #1; 2; 3  dispersion parameter


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

####ad-hoc noise
noise.mu <- 10
noise.theta <- 10
noise <- matrix(rnegbin(T0*I*simu, noise.mu, noise.theta), nrow = T0*I, ncol = simu)
Reward <- Reward + noise

####dataset####
Data <- list(Mot=Mot, Feed=Feed, Reward=Reward, id=id, t=t0, true=true_para)

####plot####
# library(ggplot2)
# colMeans(Reward)
# data_frame <- data.frame(Mot[,1], Feed[,1], Reward[,1], id[,1], t0[,1])
# ggplot(data_frame, aes(x = id[,1], y = Reward[,1])) + geom_point() + labs(x="Subject",y="Step count")
# ggplot(data_frame, aes(x = t0[,1], y = Reward[,1])) + geom_point() + labs(x="Study day",y="Step count")

####Analysis####
Estresults <- function(Data, transformation) {
  PointEst <- function(Mot, Feed, t, Reward, id, theta.temp) {
    switch(transformation,
           Poi=coef(geeglm(Reward ~ t + factor(Mot) + factor(Feed), family = poisson, id = id, corstr = "independence"))[-1],
           Quasipoi=coef(geem(Reward ~ t + factor(Mot) + factor(Feed), family = quasipoisson, id = id, corstr = "independence"))[-1],
           NB=coef(geem(Reward ~ t + factor(Mot) + factor(Feed), family = negative.binomial(theta=theta.temp), id = id, corstr = "independence"))[-1],
           glmPoi=coef(glmer(Reward ~ t + factor(Mot) + factor(Feed) + (1|id), family = poisson))$id[1,-1],
           glmQuasipoi=coef(glmer(Reward ~ t + factor(Mot) + factor(Feed) + (1|id), family = quasipoisson))$id[1,-1],
           glmNB=coef(glmer.nb(Reward ~ t + factor(Mot) + factor(Feed) + (1|id)))$id[1,-1],
           log=coef(lmer(log(Reward + 1) ~ t + factor(Mot) + factor(Feed) + (1|id)))$id[1,-1],
           SR=coef(lmer(sqrt(Reward) ~ t + factor(Mot) + factor(Feed) + (1|id)))$id[1,-1],
           bcox=coef(lmer(((Reward^0.5 - 1)/0.5) ~ t + factor(Mot) + factor(Feed) + (1|id)))$id[1,-1])
  }
  Est <- matrix(nrow = 8, ncol = simu)
  Bias <- matrix(nrow = 8, ncol = simu)
  Relative.bias <- matrix(nrow = 8, ncol = simu)
  
  for (s in (1:simu)) {
    theta.temp <- theta.ml(Data$Reward[,s], fitted(glm.nb(Data$Reward[,s] ~ Data$t[,s] + factor(Data$Mot[,s]) + factor(Data$Feed[,s]))))
    Est[,s]=unlist(PointEst(Data$Mot[,s], Data$Feed[,s], Data$t[,s], Data$Reward[,s], Data$id[,s], theta.temp))
    Bias[,s]=Est[,s] - Data$true
    Relative.bias[,s]=Bias[,s]/Data$true
  }
  
  return(
  data.frame(
  true=Data$true,
  Meanest=apply(Est, 1, mean),
  EmpiricalSE=apply(Est, 1, sd),
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
    Quasipoi=Estresults(Data, transformation="Quasipoi"),
    #glmPoi=Estresults(Data, transformation="glmPoi"), 
    #glmQuasipoi=Estresults(Data, transformation="glmQuasipoi"),
    #glmNB=Estresults(Data, transformation="glmNB"),
    #SR=Estresults(Data, transformation="SR"),
    log=Estresults(Data, transformation="log")
    #bcox=Estresults(Data, transformation="bcox")
   )
}

####Comparison results
results <- GetAnalyses(Data)

