library(MASS)
# Set the working directory (so you know where all hte files will end up!)
setwd("C:/Bob/Work/Manuscripts/Johan Standardization/log Transformation problem/Resubmission")


###
### Functions
###

# Estimate bias & RMSE for data
# transformation: which model to use
# Data includes true value, as item in the list
EstStats=function(Data, transformation, Add=0) {
  EstBias=function(y, x, True, add=0) {
    require(MASS)
    switch(transformation,
           log=coef(lm(log(y+add)~factor(x)-1))-log(True),
           sqrt=log(coef(lm(sqrt(y+add)~factor(x)-1))^2)-log(True),
           negbin=coef(glm.nb(y~factor(x)-1))-log(True),
           qpois=coef(glm(y~factor(x)-1, family=quasipoisson()))-log(True),
           coef(lm((y+add)~factor(x)-1))-log(True+add))
  }
  Bias=apply(Data$y, 2, EstBias, x=Data$Mean, True=Data$True, add=Add)
  return(data.frame(
              True=Data$True,
              MeanBias=apply(Bias, 1, mean), 
              MeanRMSE=apply(Bias, 1, function(x) sqrt(mean(x^2)))
              )
         )
}

# Function to fit all models
GetAnalyses=function(Data) {
  list(
    Data=Data,
    NB=EstStats(Data, transformation="negbin"),
    QPois=EstStats(Data, transformation="qpois"),
    Sqrt=EstStats(Data, transformation="sqrt"),
    Log0001=EstStats(Data, transformation="log", Add=0.0001),
    Log01=EstStats(Data, transformation="log", Add=0.1),
    Log05=EstStats(Data, transformation="log", Add=0.5),
    Log1=EstStats(Data, transformation="log", Add=1)
  )
}

# Plot the results from GetAnalyses
PlotStats=function(an, wh="MeanBias", ...) {
  plot(an$NB$True,an$NB[[wh]], type="l", las=1, ann=F, lty=1,...)
   lines(an$Sqrt$True,an$Sqrt[[wh]], lty=2, col=1)
#   lines(an$QPois$True,an$QPois[[wh]], lty=1, col=4, lwd=1)

   lines(an$Log1$True,an$Log1[[wh]], lty=1, col=2)
   lines(an$Log05$True,an$Log05[[wh]], lty="84", col=2)
   lines(an$Log01$True,an$Log01[[wh]], lty="34", col=2)
   lines(an$Log0001$True,an$Log0001[[wh]], lty="14", col=2)
}

############################################

###
### Simulate the data
###

True=1:20
Mean=rep(True,each=100)
NRep=500

# True: True=true mean for each level
# Mean: Mean for each data point (change to factor later!)
# y: Simulated data.  Each column is a data set.

DataHalf=list(True=True, Mean=Mean, y=replicate(NRep,rnbinom(length(Mean), 0.5, mu=Mean)))
DataOne=list(True=True, Mean=Mean, y=replicate(NRep,rnbinom(length(Mean), 1, mu=Mean)))
DataTwo=list(True=True, Mean=Mean, y=replicate(NRep,rnbinom(length(Mean), 2, mu=Mean)))
DataFive=list(True=True, Mean=Mean, y=replicate(NRep,rnbinom(length(Mean), 5, mu=Mean)))
DataTen=list(True=True, Mean=Mean, y=replicate(NRep,rnbinom(length(Mean), 10, mu=Mean)))
DataHundred=list(True=True, Mean=Mean, y=replicate(NRep,rnbinom(length(Mean), 100, mu=Mean)))

DataNames=c("DataHalf", "DataOne", "DataTwo", "DataFive", "DataTen", "DataHundred")
DataR=c(0.5, 1,2,5,10,100)

# Run all the analyses
Analyses=sapply(DataNames, function(name) {  GetAnalyses(eval(as.name(name)))  }, simplify=F)

###
### Plot the results
###

# Plot the proportion of zeroes in the data
GetP=function(data) {
  PZero=aggregate(data$y, list(Mean=data$Mean), function(lst) mean(c(lst==0)))
    PZero$p=apply(PZero[,grep("V",names(PZero))],1,mean)
  return(data.frame(Mean=PZero$Mean, P=PZero$p))
}

DataHalf.PZero=GetP(DataHalf)
DataOne.PZero=GetP(DataOne)
DataTwo.PZero=GetP(DataTwo)
DataFive.PZero=GetP(DataFive)
DataTen.PZero=GetP(DataTen)
DataHundred.PZero=GetP(DataHundred)

png("Fig1.png", height=360)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(DataHalf.PZero$Mean, DataHalf.PZero$P, type="l", ylim=c(0,0.6), tcl=0.5,
      xlab="Mean", ylab="Proportion of zeroes")
  lines(DataOne.PZero$Mean, DataOne.PZero$P, lty="94")
  lines(DataTwo.PZero$Mean, DataTwo.PZero$P, lty="74")
  lines(DataFive.PZero$Mean, DataFive.PZero$P, lty="54")
  lines(DataTen.PZero$Mean, DataTen.PZero$P, lty="34")
  lines(DataHundred.PZero$Mean, DataHundred.PZero$P, lty="14")
legend(12,0.58, as.expression(lapply(DataR, function(m) bquote(theta==.(m)))), lty=c("solid", "94", "74", "54", "34", "14"))

dev.off()

# Plot bias
png("Fig2.png", width=540, height=420)
par(mfrow=c(2,3), mar=c(3,3,2,1), oma=c(2,2,0,0))
 sapply(1:length(DataNames), function(ind, An, number, ...) {
    PlotStats(Analyses[[ind]], tcl=0.5, ...)
#    mtext(paste("r = ",number[ind],sep=""), 3,line=0.5)
    mtext(substitute(theta == num, list(num=number[ind])), 3,line=0.5)
}, An=Analyses, number=DataR, ylim=c(-3,0.5))
  legend(5.2,-1, c("Neg Bin", "Sqrt", "log, +1", "log, +0.5", "log, +0.1", "log, +0.001"),
          lty=c("solid","44", "solid","84", "34", "14"), col=c(1,1,rep(2,4)), ncol=1)
   mtext("Bias",2, outer=T)
   mtext("True Mean",1, outer=T)
dev.off()

# Plot RMSE
png("Fig3.png", width=540, height=420)
par(mfrow=c(2,3), mar=c(3,3,2,1), oma=c(2,2,0,0))
 sapply(1:length(DataNames), function(ind, An, number, ...) {
    PlotStats(Analyses[[ind]],wh="MeanRMSE", tcl=0.5, ...)
    mtext(paste("r = ",number[ind],sep=""), 3,line=0.5)
}, An=Analyses, number=DataR, ylim=c(0,1.5), yaxs="i")
  legend(5.2,1, c("Neg Bin", "Sqrt", "log, +1", "log, +0.5", "log, +0.1", "log, +0.001"),
          lty=c("solid","44","solid","84", "34", "14"), col=c(1,1,rep(2,4)), ncol=1)

   mtext("RMSE",2, outer=T)
   mtext("True Mean",1, outer=T)
dev.off()


###
### Compare negative binomial and quasipoisson models
###

diffBias= lapply(Analyses, function(lst) return(lst[["NB"]]$MeanBias-lst[["QPois"]]$MeanBias))
  range(unlist(diffBias))

diffRMSE= lapply(Analyses, function(lst) return(lst[["NB"]]$MeanRMSE-lst[["QPois"]]$MeanRMSE))
  range(unlist(diffRMSE))

plot(Analyses[["DataHalf"]][["NB"]]$True, Analyses[["DataHalf"]][["NB"]]$MeanBias-Analyses[["DataHalf"]][["QPois"]]$MeanBias, type="n",
         ylim=c(-1.0E-12,1.0E-16), xlab="True Mean", ylab="Difference in bias")
  sapply(1:length(Analyses), function(ind, List, LTY) {
     lst=List[[ind]]
     lines(lst[["NB"]]$True, lst[["NB"]]$MeanBias-lst[["QPois"]]$MeanBias, lty=LTY[ind])
  }, List=Analyses, LTY=c("solid", "94", "74", "54", "34", "14")
)