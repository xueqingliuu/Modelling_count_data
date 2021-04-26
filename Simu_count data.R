library(purrr)
library(geepack)
library(MASS)
library(lme4)
library(geeM)


####parameters####
simu <- 500
I = 100 #number of participants
# I = rdunif(1, 10000, 2000)
T0 = 45 #number of follow-ups

sigma.a <- 0.1 # 0.5 
alpha <- rnorm(I, 0, sigma.a^2) #subject-specific intercept
beta.t = 0.495
beta.M = c(0, 5.53, 6.03, 6.02) 
beta.F = c(0, 3, 3.21, 4.3, -2.5) #how to specify?
true_para <- c(beta.t, beta.M[-1], beta.F[-1])
#a <- 0.25 #1;25  overdispersion parameter

####Output matrix####
Poi.est <- matrix(nrow = 8, ncol = simu)
NB.est <- matrix(nrow = 8, ncol = simu)
log.est <- matrix(nrow = 8, ncol = simu)
SR.est <- matrix(nrow = 8, ncol = simu)
#bcox.est <- matrix(nrow = 8, ncol = simu)

Poi.bias <- matrix(nrow = 8, ncol = simu)
NB.bias <- matrix(nrow = 8, ncol = simu)
log.bias <- matrix(nrow = 8, ncol = simu)
SR.bias <- matrix(nrow = 8, ncol = simu)
#bcox.bias <- matrix(nrow = 8, ncol = simu)

Poi.relbias <- matrix(nrow = 8, ncol = simu)
NB.relbias<- matrix(nrow = 8, ncol = simu)
log.relbias <- matrix(nrow = 8, ncol = simu)
SR.relbias <- matrix(nrow = 8, ncol = simu)
#bcox.relbias <- matrix(nrow = 8, ncol = simu)


####Simulation####
for (s in 1:simu) {
  #Initial
  Reward <- rep(NA, T0*I)
  Mot <- rep(NA, T0*I)
  Feed <- rep(NA, T0*I)
  t0 <- rep(seq(1, T0, by=1), I)
  id <- rep(seq(1, I, by=1), each = 45)
  
  ####Simulate data####
  for (i in 1:I) {
    alpha.i <- alpha[i]
    for (t in 1:T0) {
      ####independent variable####
      Mot.ti <- rdunif(1, 4, 1) #Motivational message of user i at time t
      Feed.ti <- rdunif(1, 5, 1) #Feedback message of user i at time t
      Mot[(i-1)*T0 + t] <- Mot.ti
      Feed[(i-1)*T0 + t] <- Feed.ti
      ####dependent variable#### (Possion)
      lamb.ti <- exp(alpha.i + beta.t*t + beta.M[Mot.ti] + beta.F[Feed.ti])
      #theta.ti <- rgamma(1, shape = a, rate = 1/a) #overdispersion
      theta.ti <- 1 #without overdispersion
      Reward[(i-1)*T0 + t] <- rpois(1, lamb.ti*theta.ti)
    }
  }
  
  Mot <- as.factor(Mot)
  Feed <- as.factor(Feed)
  data <- data.frame(Mot, Feed, id, t, Reward)
  ####Possion Regression####
  Poi.geeglm <- geeglm(Reward ~ t0 + Mot + Feed, family = poisson, id = id, corstr = "independence")
  Poi.est[,s] <- coef(Poi.geeglm)[-1] #estimates
  Poi.bias[,s] <- Poi.est[,s] - true_para #bias
  Poi.relbias[,s] <- Poi.bias[,s]/true_para #relative bias

  ####Quasi-Possion Regression####
  #QuasiPoi.glm <- glm(Reward ~ t0 + Mot + Feed, family = quasipoisson)
  #QuasiPoi.geeglm <- geeglm(Reward ~ t0 + Mot + Feed, family = quasipoisson, id = id, corstr = "independence")
  
  
  ####Negative Binomial Regression####
  NB.geeglm <- geem(Reward ~ t0 + Mot + Feed, family = negative.binomial(theta=2), id = id, corstr = "independence")
  NB.est[,s] <- coef(NB.geeglm)[-1]
  NB.bias[,s] <- NB.est[,s] - true_para
  NB.relbias[,s] <- NB.bias[,s]/true_para

  #NB.geeglm <- geeglm(Reward ~ t0 + Mot + Feed, family = negative.binomial(theta=2), id = id, corstr = "independence")
  
  ####Log Transformation####
  log.lmer <- lmer(log(Reward) ~ t0 + Mot + Feed + (1|id))
  log.est[,s] <- coef(log.lmer)$id[1,-1]
  log.bias[,s] <- log.est[,s] - true_para
  log.relbias[,s] <- log.bias[,s]/true_para
  
  ####Square Root Transformation#### (how to compare it with others?)
  SR.lmer <- lmer(sqrt(Reward) ~ t0 + Mot + Feed + (1|id))
  SR.est[,s] <- coef(SR.lmer)$id[1,-1]
  SR.bias[,s] <- SR.est[,s] - true_para
  SR.relbias[,s] <- SR.bias[,s]/true_para
  
  ####Box-Cox Transformation####
  # lam = 0.5
  # bcox.lmer <- lmer(((Reward^lam-1)/lam) ~ t0 + Mot + Feed + (1|id))
  # bcox.est[,s] <- coef(bcox.lmer$id[1,-1])
  # bcox.bias[,s] <- bcox.est[,s] - true_para
  # bcox.relbias[,s] <- bcox.est[,s]/true_para
  
}


####Performance measure results####
#Poisson
apply(Poi.est, 1, mean) #mean estimates
apply(Poi.est, 1, var) #variance
apply(Poi.bias, 1, mean) #mean mse
apply(Poi.bias, 1, function(x) mean(x^2)) #MSE
apply(Poi.relbias, 1, function(x) mean) # mean relative bias
#NB
apply(NB.est, 1, mean) #mean estimates
apply(NB.est, 1, var) #variance
apply(NB.bias, 1, mean) #mean mse
apply(NB.bias, 1, function(x) mean(x^2)) #MSE
apply(NB.relbias, 1, function(x) mean) # mean relative bias
#log
apply(log.est, 1, mean) #mean estimates
apply(log.est, 1, var) #variance
apply(log.bias, 1, mean) #mean mse
apply(log.bias, 1, function(x) mean(x^2)) #MSE
apply(log.relbias, 1, function(x) mean) # mean relative bias
#SR
apply(SR.est, 1, mean) #mean estimates
apply(SR.est, 1, var) #variance
apply(SR.bias, 1, mean) #mean mse
apply(SR.bias, 1, function(x) mean(x^2)) #MSE
apply(SR.relbias, 1, function(x) mean) # mean relative bias
#boxcox
# apply(bcox.est, 1, mean) #mean estimates
# apply(bcox.est, 1, var) #variance
# apply(bcox.bias, 1, mean) #mean mse
# apply(bcox.bias, 1, function(x) mean(x^2)) #MSE
# apply(bcox.relbias, 1, function(x) mean) # mean relative bias
