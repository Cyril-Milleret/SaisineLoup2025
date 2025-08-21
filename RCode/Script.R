rm(list=ls())
setwd("C:/Personal_Cloud/OneDrive/Work/CNRS/ProjectionHarvest/Esco2025/GitHub/SaisineLoup2025/RCode")

#loacd packages
library(nimble)
library(tidyverse)
library(xtable)
library(basicMCMCplots)
library(coda)
library(readxl)
library(dplyr)
library(raster)
#Load a few functions necessary for plotting and rearanging the data
source("ProcessCodaOutput.R")
#modified version of the basicMCMCplots::chainsPlot() to print the output as picture
source("chainsPlot1.R")

#working directory where figures will be stored 
wdFigure <- "C:/Personal_Cloud/OneDrive/Work/CNRS/ProjectionHarvest/Esco2025/GitHub/SaisineLoup2025/Figures"

#set chain characteristics
n.iter = 200000
n.burnin = 20000
nchains = 4
# the data 
harvest <- c(0,0,0,0,0,1,0,0,2,1,2,0,0,1,0,4,4,6,18,36,34,42,51,98,105,103,169,207,204)
CMR <- c(17.1,35.4,47.7,25.1,62.6,47.9,81.7,110.5,102.7,135.9,132.6,101.7,130.3,141.4,141.5,175.5,210.3,174.5,353.6,280.2,376.7,561.2,571.9,682.4,645.7,783.8,1081,1002.5,1013)
ans <- 1996:(1996 + length(harvest) - 1)
thedata <- cbind(ans,round(CMR), harvest)
colnames(thedata) <- c("year", "N", "H")
thedata <- as.data.frame(thedata)
nyears <- nrow(thedata)

#==== 1. Growth rate estimation with legal culling (Observed lambda) ==== 
#====   1.1 Exponential model with constant lambda ==== 
#====     1.1.1 model definition ==== 
modelExpNoH <- nimbleCode({
  # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/(sigmaProc^2)
  lambda ~ dunif(0, 2)
  N[1] ~ dgamma(1.0E-6, 1.0E-6)
  
  # Process model
  for (t in 2:(nyears)) {
    Nproc[t] <- log(max(1, lambda*(N[t-1])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
  }
  
})
#====     1.1.2 Prepare nimble objects and fit model ====
nimData <- list(Nobs = thedata$N,
                H = thedata$H
)

nimConstants <- list(nyears = nyears)
nimParams <- c('lambda', "sigmaProc","errorObs")
set.seed(10)
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears,10,15),
                  N = nimData$Nobs
)


## Create and compile the NIMBLE model
model <- nimbleModel( code = modelExpNoH,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                          enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
exp_resNoHSamples <- runMCMC( mcmc = cMCMC,
                              nburnin = n.burnin,
                              niter = n.iter,
                              nchains = nchains ,
                              samplesAsCodaMCMC = TRUE,
                              WAIC = T)
#get sim.lists
exp_resNoH <- ProcessCodaOutput(exp_resNoHSamples$samples)

#====     1.1.3 Plot results ====
.expressions <- c("lambda",
                  "sigma_proc",
                  "sigma_obs",
                  "shapeObs")

# "sigmaObs","errorObs","shapeObs"
legend_expressions <- parse(text = .expressions)
#check convergence 
jpeg(filename =paste(wdFigure,"/ConvergenceObservedConstant.jpg",sep=""), width = 1000,
     height = 1500,pointsize = 28,quality =500)
chainsPlot1(samplesList=exp_resNoHSamples$samples, 
            var = c('lambda', "sigmaProc"),#,
            titleName = c("lambda",
                          "sigma_proc"
                          ))
dev.off()


#====   1.2 Exponential model with period effect  ====
#====     1.2.1 model definition ==== 
modelPeriodNoH <- nimbleCode({
  # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda ~ dunif(0, 2)
  N[1] ~ dgamma(1.0E-6, 1.0E-6)
  beta ~ dnorm(0, 1.0E-6)
  
  
  # Process model
  for (t in 2:(nyears)) {
    Nproc[t] <- log(max(1, (lambda + beta*X[t-1])*(N[t-1])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
  }
})
#====     1.2.2 Prepare nimble objects and fit model  ====
# here we define two periods 
nimData <- list(Nobs = thedata$N,
                H=thedata$H,
                X = c(rep(0, nyears-6), rep(1, 5), NA)
)

## check period definition
thedata$year[nimData$X%in%0]
thedata$year[nimData$X%in%1]

nimConstants <- list(nyears = nyears)

nimParams <- c('lambda', "sigmaProc", "beta","errorObs")
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears,10,15),
                  beta =  runif(1,-0.1,0.1),
                  N = nimData$Nobs
)


## Create and compile the NIMBLE model
model <- nimbleModel( code = modelPeriodNoH,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                          enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
period_mcmcNoHsamples <- runMCMC( mcmc = cMCMC,
                                 nburnin = n.burnin,
                                 niter = n.iter,
                                 nchains = nchains ,
                                 samplesAsCodaMCMC = TRUE,
                                 WAIC = T)
#get sim.lists
period_mcmcNoH <- ProcessCodaOutput(period_mcmcNoHsamples$samples)

#====     1.2.3 Plot results ====
#check convergence 
jpeg(filename =paste(wdFigure,"/ConvergenceObservedPeriod.jpg",sep=""), width = 1000,
     height = 1500,pointsize = 28,quality =500)
chainsPlot1(samplesList=period_mcmcNoHsamples$samples, 
            var = c('lambda', "sigmaProc","beta"),#,
            titleName = c("lambda",
                          "sigma_proc",
                          "beta"
            ))
# )#,file =paste(wdFigure,"/Convergence.pdf",sep=""))
dev.off()



post_periodNOH <- period_mcmcNoH$sims.list %>%
  as_tibble() %>%
  mutate(lambda1 =lambda+beta )%>%
  #select(lambda1,lambda)%>%
  #  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  #  filter(str_detect(parameter, "lambda")) %>%
  ggplot() +
  geom_density(aes(x = lambda),color = "green4") + 
  geom_density(aes(x = lambda1),color ="purple") + 
  # geom_density(aes(x = lambda)) + 
  geom_vline(xintercept = 1, lty = "dashed", color = "red") +
  labs(x = "Taux de croissance", title = "freinage")+
  geom_text(x=1.3, y=2, label=paste("2nd Period ",round(period_mcmcNoH$mean$lambda+
                                                          period_mcmcNoH$mean$beta,digits=2), " (95%: ",
                                    round(quantile(period_mcmcNoH$sims.list$lambda+
                                                     period_mcmcNoH$sims.list$beta
                                                   ,probs = 0.025),digits=2),
                                    "-", 
                                    round(quantile(period_mcmcNoH$sims.list$lambda+
                                                     period_mcmcNoH$sims.list$beta,probs = 0.975),digits=2),
                                    ")",sep=""),color ="purple")+
  geom_text(x=1.3, y=5, label=paste("1st Period ",round(period_mcmcNoH$mean$lambda,digits=2), " (95%: ",
                                    round(quantile(period_mcmcNoH$sims.list$lambda
                                                   ,probs = 0.025),digits=2),
                                    "-", 
                                    round(quantile(period_mcmcNoH$sims.list$lambda,
                                                   probs = 0.975),digits=2),
                                    ")",sep=""),color ="green4")

# post_periodNOH




#====   1.3 Exponential model with random effect  ====
#====     1.3.1 model definition ==== 
modelRDENoH <- nimbleCode({
  # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda0 ~ dnorm(0, 0.0001) #random growth rate
  lambdaSigma ~ dunif(0,10)
  lambdaTau <- 1/(lambdaSigma^2)
  
  N[1] ~ dgamma(1.0E-6, 1.0E-6)
  
  # Process model
  for (t in 2:(nyears)) {
    lambdat[t-1] ~ dnorm(lambda0, tau= lambdaTau) # random growth rate
    Nproc[t] <- log(max(1, exp(lambdat[t-1])*(N[t-1])))#
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
  }
  
})

#====     1.3.2 Prepare nimble objects and fit model  ====
nimParams <- c("lambda0","lambdat","lambdaSigma",
               "lambdaTau", "sigmaProc", "sigmaObs",
               "N", "tauProc")

#adjust initial values 
set.seed(10)
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda0 = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears,10,15),
                  lambdaSigma= runif(1,0.1,0.2),
                  lambdat = runif(nimConstants$nyears-1,0.1,0.2),
                  lambdaPred = runif(length(nimData$harvest),0.1,0.2),
                  N = nimData$Nobs
                  
)

nimData <- list(Nobs = thedata$N,
                H=thedata$H)


nimConstants <- list(nyears = nyears)

## Create and compile the NIMBLE model
model <- nimbleModel( code = modelRDENoH,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                          enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(6)
rde_mcmcNoHsamples <- runMCMC( mcmc = cMCMC,
                                  nburnin = n.burnin,
                                  niter = n.iter,
                                  nchains = nchains ,
                                  samplesAsCodaMCMC = TRUE,
                                  WAIC = T)
#get sim.lists
rde_mcmcNoH <- ProcessCodaOutput(rde_mcmcNoHsamples$samples)

#====     1.3.3 Plot results ====
#check convergence 
jpeg(filename =paste(wdFigure,"/ConvergenceObservedRandom.jpg",sep=""), width = 1000,
     height = 1500,pointsize = 28,quality =500)
chainsPlot1(samplesList=rde_mcmcNoHsamples$samples, 
            var = c('lambda0', "sigmaProc","beta"),#,
            titleName = c("lambda0",
                          "sigma_proc"
            ))
# )#,file =paste(wdFigure,"/Convergence.pdf",sep=""))
dev.off()
#====   1.4 Compare models WAIC ==== 
round(exp_resNoHSamples$WAIC$WAIC,digits=2)
round(period_mcmcNoHsamples$WAIC$WAIC,digits=2)
round(rde_mcmcNoHsamples$WAIC$WAIC,digits=2)
#====   1.5 Plot lambda all three models together ==== 
png(filename =  paste(wdFigure,"/lambdaObserved.png",sep=""), pointsize = 20,
    width = 1380,height = 800)
par(mfrow=c(1,3))

#constant
dens <- density(exp_resNoH$sims.list$lambda)
plot(dens, xlim=c(0.8,1.5),
     xlab = expression(lambda[obs]), ylab="Densité", main="Constant", 
     col=grey(0.2),cex.axis=1.2,cex.lab=1.2)
polygon(dens, col = grey(0.2,alpha=0.5))
abline(v=1,col="red",lty=2)

text(exp_resNoH$mean$lambda+0.18,11,paste(round(exp_resNoH$mean$lambda,digits=2), " (",
           round(quantile(exp_resNoH$sims.list$lambda,probs = 0.025),digits=2),
           "-", 
           round(quantile(exp_resNoH$sims.list$lambda,probs = 0.975),digits=2),
           ")",sep=""),cex=1.5)

dens <- density(period_mcmcNoH$sims.list$lambda)

#period effect
plot(dens, xlim=c(0.8,1.5),
     xlab = expression(lambda[obs]), ylab="Densité", main="Effet période",cex.axis=1.2,cex.lab=1.2)
polygon(dens, col = grey(0.2,alpha=0.5))
abline(v=1,col="red",lty=2)
lambda1 <- period_mcmcNoH$sims.list$lambda+period_mcmcNoH$sims.list$beta
dens <- density(lambda1)
lines(dens,
      col=grey(0.8))
polygon(dens, col = grey(0.8,alpha=0.5))
text(period_mcmcNoH$mean$lambda+0.18,10,paste("P1<2019 \n",format(round(period_mcmcNoH$mean$lambda,digits=2),nsmall = 2), " (",
                                              format(round(quantile(period_mcmcNoH$sims.list$lambda,probs = 0.025),digits=2),nsmall = 2),
                                          "-", 
                                              format(round(quantile(period_mcmcNoH$sims.list$lambda,probs = 0.975),digits=2),nsmall = 2),
                                          ")",sep=""),cex=1.5)

text(period_mcmcNoH$mean$lambda+0.22,4,paste("P2", "≥","2019","\n",format(round(mean(lambda1),digits=2),nsmall = 2), " (",
                                             format(round(quantile(lambda1,probs = 0.025),digits=2),nsmall = 2),
                                             "-", 
                                             format(round(quantile(lambda1,probs = 0.975),digits=2),nsmall = 2),
                                             ")",sep=""),cex=1.5)
lines(x=c(1.16,period_mcmcNoH$mean$lambda+0.09),y=rep(4,2))

#random effect
dens <- density(exp(rde_mcmcNoH$sims.list$lambda0))
plot(dens, xlim=c(0.8,1.5),
     xlab = expression(lambda[obs]), ylab="Densité", main="Effet aléatoire",cex.axis=1.2,cex.lab=1.2)
abline(v=1,col="red",lty=2)
polygon(dens, col = grey(0.2,alpha=0.5))
text(exp(rde_mcmcNoH$mean$lambda0)+0.18,11,paste(round(exp(rde_mcmcNoH$mean$lambda0),digits=2), " (",
                                           round(quantile(exp(rde_mcmcNoH$sims.list$lambda0),probs = 0.025),digits=2),
                                           "-", 
                                           round(quantile(exp(rde_mcmcNoH$sims.list$lambda0),probs = 0.975),digits=2),
                                           ")",sep=""),cex=1.5)
dev.off()

#==== 2. Growth rate estimation and projection accounting for legal culling (potential lambda) ==== 
#Define range of legal culling used to project  N at t+1
harvestPred <- 0:500
#add to nimble object
nimData$harvest <- harvestPred
nimConstants$NharvestPred <- length(harvestPred)
#predict until 2035/36
NyearsPredic <- 2036 - thedata$year[nyears]
NyearsPred <- nyears + NyearsPredic
nimConstants$NyearsPred <- NyearsPred

#====   2.1 Exponential model with H  ====
modelExp <- nimbleCode({
 # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda ~ dunif(0, 2)
  N[1] ~ dgamma(1.0E-6, 1.0E-6)

  # Process model
  for (t in 2:(nyears)) {
    Nproc[t] <- log(max(1, lambda*(N[t-1]- H[t-1])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
}

  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
}
  
  # Projected population until 2035
  for (t in (nyears + 1):NyearsPred) {
    #here we fix the rate to 19%
    Nkilled[t] <- N[t-1] *0.19
    Nproc[t] <- log(max(1, lambda*(N[t-1]-Nkilled[t])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }

##################################
# Harvested population
##################################
# here we make prediction for t +1 for different number of individuals killed
  for (i in 1:NharvestPred) {
	Nprocharv[i] <- log(max(1, lambda*(N[nyears] - harvest[i]))) # Because of index
	Ntrueharv[i] ~ dlnorm(Nprocharv[i], tauProc)
	growth[i] <- Ntrueharv[i]/N[nyears]
}
  
})

#====   2.2 Exponential model with H and period effect on lambda====

modelPeriod <- nimbleCode({
  # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda ~ dunif(0, 2)
  beta ~ dnorm(0, 1.0E-6)
  N[1] ~ dgamma(1.0E-6, 1.0E-6)
  
  # Process model
  for (t in 2:(nyears)) {
    Nproc[t] <- log(max(1, (lambda + beta*X[t-1])*(N[t-1]- H[t-1])))#
    N[t] ~ dlnorm(Nproc[t], tauProc)
}
  
  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
}
  
  # Projected population
  for (t in (nyears + 1):NyearsPred) {
    #here we fix the rate to 19%
    Nkilled[t] <- N[t-1] * 0.19
    Nproc[t] <- log(max(1, (lambda + beta)*(N[t-1]-Nkilled[t])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
##################################
# Harvested population
##################################
  for (i in 1:NharvestPred) {
	procNharv[i] <- log(max(1, (lambda+beta)*(N[nyears] - harvest[i]))) # Because of index
	Ntrueharv[i] ~ dlnorm(procNharv[i], tauProc)
	growth[i] <- Ntrueharv[i]/N[nyears]
  }
  
}
)

#====   2.3 Exponential model with H and random effect on lambda ====

modelExpRde <- nimbleCode({
 # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/(sigmaProc^2)
  lambda0 ~ dnorm(0, 0.0001) 
  lambdaSigma ~ dunif(0,10)
  lambdaTau <- 1/(lambdaSigma^2)
  
  N[1] ~ dgamma(1.0E-6, 1.0E-6)

  # Process model
  for (t in 2:(nyears)) {
    lambdat[t-1] ~ dnorm(lambda0, tau= lambdaTau) # random growth rate
    Nproc[t] <- log(max(1, exp(lambdat[t-1])*(N[t-1]-H[t-1])))#H
    N[t] ~ dlnorm(Nproc[t], tauProc)
}

  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
}
  
  # Projected population until 2035
  for (t in (nyears + 1):(NyearsPred)) {
    #here we fix the rate to 19%
    Nkilled[t] <- N[t-1] *0.19
    lambdat[t-1] ~ dnorm(lambda0, tau= lambdaTau) # random growth rate

    Nproc[t] <- log(max(1, exp(lambdat[t-1])*(N[t-1]-Nkilled[t])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
  
# Harvested population
# here we make prediction for t +1 for different number of individuals killed
  for (i in 1:NharvestPred) {
  lambdaPred[i] ~ dnorm(lambda0, tau= lambdaTau) # random growth rate
	Nprocharv[i] <- log(max(1, exp(lambdaPred[i])*(N[nyears] - harvest[i]))) # Because of index
	Ntrueharv[i] ~ dlnorm(Nprocharv[i], tauProc)
	growth[i] <- Ntrueharv[i]/N[nyears]
}
  
})


#====   2.4 Prepare nimble objects and fit nimble model plot results ==== 
#====     2.4.1 Exponential model  ====
#====       2.4.1.1 Prepare nimble objects and fit model====
nimParams <- c("lambda", "sigmaProc", "sigmaObs", "N", "tauProc","Ntrueharv", "growth","errorObs")

#adjust initial values 
set.seed(10)
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears,10,15),
                  beta =  runif(1,-0.1,0.1),
                  N = nimData$Nobs
)
nimInits$N <- c(nimInits$N, rep(nimInits$N[nyears],NyearsPredic))
nimInits$lambdaNobs <- c(nimInits$lambdaNobs, rep(nimInits$lambdaNobs[nyears], NyearsPredic))
nimInits$Ntrueharv <- c(rep(nimInits$lambdaNobs[nyears], length(nimData$harvest)))
nimInits$Nkilled <- c(nimInits$lambdaNobs,rep(nimInits$lambdaNobs[nyears], NyearsPredic))


## Create and compile the NIMBLE model
model <- nimbleModel( code = modelExp,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                          enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
exp_mcmcsamples <- runMCMC( mcmc = cMCMC,
                            nburnin = n.burnin,
                            niter = n.iter,
                            nchains = nchains ,
                            samplesAsCodaMCMC = TRUE,
                            WAIC = T)
#get sim.lists
exp_mcmc <- ProcessCodaOutput(exp_mcmcsamples$samples,params.omit =c("Ntrueharv","sigmaObs","growth","N"),DIC=F )


exp_mcmc$Rhat
#====       2.4.1.2 Plot results====
#check convergence 
jpeg(filename =paste(wdFigure,"/ConvergencePotentialConstant.jpg",sep=""), width = 1000,
     height = 1500,pointsize = 28,quality =500)
chainsPlot1(samplesList=exp_mcmcsamples$samples, 
            var = c('lambda', "sigmaProc","errorObs"),#,
            titleName = c("lambda",
                          "sigma_proc","errorObs")
            )
# )#,file =paste(wdFigure,"/Convergence.pdf",sep=""))
dev.off()


#Projection 2035
year <- c(1996:2036)
year1 <- year-1
Winter <- paste(substr(year1, 3,4),"/",substr(year, 3,4),sep="")


dat <- exp_mcmc$sims.list$N %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(medianN = median(value),
            lci = quantile(value, probs = 2.5/100),
            uci = quantile(value, probs = 97.5/100),
            lci50 = quantile(value, probs = 25/100),
            uci50 = quantile(value, probs = 75/100),
            lci75 = quantile(value, probs = 12.5/100),
            uci75 = quantile(value, probs = 87.5/100),
            ) %>%
  mutate(an = parse_number(parameter)) %>%
  arrange(an) 


#
proj_exp2035Line5075 <- 
  ggplot() + 
  geom_ribbon(data= dat, aes(x = an, y = medianN, ymin = lci, ymax = uci), 
              fill = grey(0.6), alpha = 0.4) + 
  geom_ribbon(data= dat, aes(x = an, y = medianN, ymin = lci75, ymax = uci75), 
              fill = grey(0.4), alpha = 0.4) +
  geom_ribbon(data= dat, aes(x = an, y = medianN, ymin = lci50, ymax = uci50), 
              fill = grey(0.1), alpha = 0.4) + 
 
  # geom_ribbon(data= dat[nyears: nrow(dat),], aes(x = an, y = medianN, ymin = lci50, ymax = uci50), 
  #             fill = grey(0.3), alpha = 0.6) + 
  #geom_line(data= dat[1:nyears+0,], aes(x = an, y = medianN), lty = "dashed", color = "red") + 
  #  geom_point(aes(x = an, y = medianN), color = "red") +
  geom_point(data = thedata %>% as_tibble, aes(x = 1:unique(nyears), y = N)) + 
  coord_cartesian(xlim = c(1, 40), ylim = c(0, 4200)) +
  labs(y = "Effectifs",
       x = "Hivers")+
  #title = "Effectifs projetes selon modèle exponentiel",
  #subtitle = "avec effectifs observes (points noirs)") + 
  scale_x_continuous(breaks = seq(1, 41, by = 4),
                     labels = Winter[seq(1, 41, by = 4)])


pdf(file =  paste(wdFigure,"/Proj2035LineExp5075.pdf",sep=""),
    width = 8,height = 8)
proj_exp2035Line5075
dev.off()

#====     2.4.2. Exponential model with period effect  ====
#====       2.4.2.1 Prepare nimble objects and fit model ====
nimParams <- c("lambda", "sigmaProc", "sigmaObs", "N",
               "tauProc","Ntrueharv", "growth","beta","errorObs")
nimData$X = c(rep(0, nyears-6), rep(1, 5), NA)

model <- nimbleModel( code = modelPeriod,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                         enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
period_mcmcsamples <- runMCMC( mcmc = cMCMC,
                                 nburnin = n.burnin,
                                 niter = n.iter,
                                 nchains = nchains ,
                                 samplesAsCodaMCMC = TRUE,
                                 WAIC = T)
#get sim.lists
period_mcmc <- ProcessCodaOutput(period_mcmcsamples$samples,params.omit =c("Ntrueharv","sigmaObs","growth","N"), DIC=F )

#====       2.4.2.2 Plot results ====
jpeg(filename =paste(wdFigure,"/ConvergencePotentialPeriod.jpg",sep=""), width = 1000,
     height = 1500,pointsize = 28,quality =500)
chainsPlot1(samplesList=period_mcmcsamples$samples, 
            var = c('lambda', "sigmaProc","errorObs","beta"),#,
            titleName = c("lambda",
                          "sigma_proc","errorObs","beta")
)
# )#,file =paste(wdFigure,"/Convergence.pdf",sep=""))
dev.off()




#projection
proj_period2035 <-  period_mcmc$sims.list$N %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  #filter(str_detect(parameter, "N\\[")) %>%
  #filter(parameter== "V") %>%
  
  group_by(parameter) %>%
  summarize(medianN = median(value),
            lci = quantile(value, probs = 2.5/100),
            uci = quantile(value, probs = 97.5/100),
            lci50 = quantile(value, probs = 25/100),
            uci50 = quantile(value, probs = 75/100),
            lci75 = quantile(value, probs = 12.5/100),
            uci75 = quantile(value, probs = 87.5/100)) %>%
  mutate(an = parse_number(parameter)) %>%
  arrange(an) %>%
  ggplot() + 
  geom_ribbon(aes(x = an, y = medianN, ymin = lci, ymax = uci),  fill = grey(0.6), alpha = 0.4) + 
  geom_ribbon(aes(x = an, y = medianN, ymin = lci50, ymax = uci50), 
              fill = grey(0.1), alpha = 0.4)+
geom_point(data = thedata %>% as_tibble, aes(x = 1:unique(nyears), y = N)) + 
  coord_cartesian(xlim = c(1, 40), ylim = c(0, 4200)) +
  labs(y = "Effectifs",
       x = "Hivers")+
  #title = "Effectifs projetes selon modèle exponentiel",
  #subtitle = "avec effectifs observes (points noirs)") + 
  scale_x_continuous(breaks = seq(1, 41, by = 4),
                     labels = Winter[seq(1, 41, by = 4)])
proj_period2035



#LINE THAT STOPS AFTER FEW YEARS OF PREDICTIONS
# pdf(file =  paste(wdFigure,"/Proj2035Period.pdf",sep=""),
#     width = 8,height = 8)
# proj_period2035
# dev.off()



#====     2.4.3. Exponential model with random effect  ====
#====       2.4.3.1 Prepare nimble objects and fit model====
nimParams <- c("lambda0","lambdat","lambdaSigma",
               "lambdaTau", "sigmaProc", "sigmaObs","errorObs", "N","Ntrueharv", "growth")

#adjust initial values 
set.seed(10)
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda0 = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears,10,15),
                  N = nimData$Nobs,
                  lambdaSigma= runif(1,0.1,0.2),
                  lambdaSigma1= runif(1,0.1,0.2),
                  lambdaTau = runif(1,0.1,0.2),
                  lambdat = runif(nimConstants$nyears-1,0.1,0.2),
                  lambdaPred = runif(length(nimData$harvest),0.1,0.2)
)
nimInits$N <- c(nimInits$N, rep(nimInits$N[nyears],NyearsPredic))
nimInits$lambdaNobs <- c(nimInits$lambdaNobs, rep(nimInits$lambdaNobs[nyears], NyearsPredic))
nimInits$lambdat <- c(nimInits$lambdat, rep(nimInits$lambdat[nyears-1], NyearsPredic))
nimInits$lambda <- c(rep(nimInits$lambdat[nyears-1], length(nimData$harvest)))

nimInits$Ntrueharv <- c(rep(nimInits$lambdaNobs[nyears], length(nimData$harvest)))
nimInits$Nkilled <- c(nimInits$lambdaNobs,rep(nimInits$lambdaNobs[nyears], NyearsPredic))

model <- nimbleModel( code = modelExpRde,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          # control = list(reflective = TRUE),
                          thin = 1,
                          enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
expRde_mcmcsamples <- runMCMC( mcmc = cMCMC,
                               nburnin = n.burnin,
                               niter = n.iter*2,
                               nchains = nchains ,
                               samplesAsCodaMCMC = TRUE,
                               WAIC = T)
#get sim.lists
expRde_mcmc <- ProcessCodaOutput(expRde_mcmcsamples$samples,
                                 params.omit =c("lambda0","lambdat","lambda","Ntrueharv","tau","growth","N"),DIC=F )
#====       2.4.3.2 Plot results ====
expRde_mcmc$Rhat#does not converge
#growth
jpeg(filename =paste(wdFigure,"/ConvergencePotentialRandom.jpg",sep=""), width = 1000,
     height = 1500,pointsize = 28,quality =500)
chainsPlot1(samplesList=expRde_mcmcsamples$samples, 
            var = c('lambda0', "sigmaProc","errorObs","lambdaSigma"),#,
            titleName = c("lambda0",
                          "sigma_proc","errorObs","lambdaSigma")
)
# )#,file =paste(wdFigure,"/Convergence.pdf",sep=""))
dev.off()

##projection
proj_expRde2035 <-  expRde_mcmc$sims.list$N %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(medianN = median(value),
            lci = quantile(value, probs = 2.5/100),
            uci = quantile(value, probs = 97.5/100),
            lci50 = quantile(value, probs = 25/100),
            uci50 = quantile(value, probs = 75/100),
            lci75 = quantile(value, probs = 12.5/100),
            uci75 = quantile(value, probs = 87.5/100)) %>%
  mutate(an = parse_number(parameter)) %>%
  arrange(an) %>%
  ggplot() + 
  geom_ribbon(aes(x = an, y = medianN, ymin = lci, ymax = uci),  fill = grey(0.6), alpha = 0.4) + 
  geom_ribbon(aes(x = an, y = medianN, ymin = lci50, ymax = uci50), 
              fill = grey(0.1), alpha = 0.4)+
  geom_point(data = thedata %>% as_tibble, aes(x = 1:unique(nyears), y = N)) + 
  coord_cartesian(xlim = c(1, 40), ylim = c(0, 4200)) +
  labs(y = "Effectifs",
       x = "Hivers")+
  scale_x_continuous(breaks = seq(1, 41, by = 4),
                     labels = Winter[seq(1, 41, by = 4)])

##save plot
# 
# pdf(file =  paste(wdFigure,"/Proj2035rde.pdf",sep=""),
#     width = 8,height = 8)
# proj_expRde2035
# dev.off()

#====     2.4.4 Plot results 3 models ====
#====       2.4.4.1 Lambda ====
pdf(file =  paste(wdFigure,"/lambdaPotentiel.pdf",sep=""),
    width = 12,height = 6)
par(mfrow=c(1,3))
#constant
dens <- density(exp_mcmc$sims.list$lambda)
plot(dens, xlim=c(0.8,1.8),
     xlab = expression(lambda), ylab="Densité", main="Constant", col=grey(0.2))
polygon(dens, col = grey(0.2,alpha=0.5))
abline(v=1,col="red",lty=2)

text(exp_mcmc$mean$lambda+0.28,8,paste(round(exp_mcmc$mean$lambda,digits=2), " (",
                                          round(quantile(exp_mcmc$sims.list$lambda,probs = 0.025),digits=2),
                                          "-", 
                                          round(quantile(exp_mcmc$sims.list$lambda,probs = 0.975),digits=2),
                                          ")",sep=""),cex=1.5)

dens <- density(period_mcmc$sims.list$lambda)
#period effect
plot(dens, xlim=c(0.8,1.8),
     xlab = expression(lambda), ylab="Densité", main="Effet période")
polygon(dens, col = grey(0.2,alpha=0.5))
abline(v=1,col="red",lty=2)
lambda1 <- period_mcmc$sims.list$lambda+period_mcmc$sims.list$beta
dens <- density(lambda1)
lines(dens,
      col=grey(0.8))
polygon(dens, col = grey(0.8,alpha=0.5))

 
text(period_mcmc$mean$lambda+0.28,8,paste("P1<2019 \n",format(round(period_mcmc$mean$lambda,digits=2),nsmall = 2), " (",
                                              format(round(quantile(period_mcmc$sims.list$lambda,probs = 0.025),digits=2),nsmall = 2),
                                              "-", 
                                              format(round(quantile(period_mcmc$sims.list$lambda,probs = 0.975),digits=2),nsmall = 2),
                                              ")",sep=""),cex=1.5)

text(1.6,3,paste("P2", "≥","2019","\n",format(round(mean(lambda1),digits=2),nsmall = 2), " (",
                                             format(round(quantile(lambda1,probs = 0.025),digits=2),nsmall = 2),
                                             "-", 
                                             format(round(quantile(lambda1,probs = 0.975),digits=2),nsmall = 2),
                                             ")",sep=""),cex=1.5)

expression("T ">="5")
dens <- density(exp(expRde_mcmc$sims.list$lambda0))
plot(dens, xlim=c(0.8,1.8),
     xlab = expression(lambda), ylab="Densité", main="Effet aléatoire")
abline(v=1,col="red",lty=2)
polygon(dens, col = grey(0.2,alpha=0.5))

text(exp(expRde_mcmc$mean$lambda0)+0.28,8,paste(round(exp(expRde_mcmc$mean$lambda0),digits=2), " (",
                                        round(quantile(exp(expRde_mcmc$sims.list$lambda0),probs = 0.025),digits=2),
                                        "-", 
                                        round(quantile(exp(expRde_mcmc$sims.list$lambda0),probs = 0.975),digits=2),
                                        ")",sep=""),cex=1.5)
dev.off()

#==== 3. Model validation and choice ====
#====   3.1 WAIC ====
#exponentiel
round(exp_mcmcsamples$WAIC$WAIC,digits=2)#
#Periode
round(period_mcmcsamples$WAIC$WAIC,digits=2)
#random effect
round(expRde_mcmcsamples$WAIC$WAIC,digits=2)

#====   3.2 predictive power (Remove and Predict last two years) ====
nimData <- list(Nobs = thedata$N,
                H=thedata$H,
	            X = c(rep(0, nyears-8), rep(1, 7), NA)
)

nimConstants <- list(nyears = nyears)

set.seed(10)
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears,10,15),
                  beta =  runif(1,-0.1,0.1),
                  N = nimData$Nobs
)

NyearsPred <- 2
# We set NA for the two last years for predictions
thedata1 <- thedata
thedata1$N[(length(thedata1$N)-NyearsPred+1): length(thedata1$N)] <- NA

nimData$Nobs <- thedata1$N
nimInits$Nobs  <- thedata$N
nimInits$Nobs[!is.na(nimData$Nobs)] <- NA
nimConstants$nyears <- nimConstants$nyears-2
nimParams <- c("N")
nimInits$N <-  thedata$N

#====     3.2.1 exp model ====
modelExp <- nimbleCode({
 # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda ~ dunif(0, 2)
  N[1] ~ dgamma(1.0E-6, 1.0E-6)

  # Process model
  for (t in 2:(nyears)) {
    Nproc[t] <- log(max(1, lambda*(N[t-1]- H[t-1])))#
    N[t] ~ dlnorm(Nproc[t], tauProc)
}

  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
}
  # Projected population until 2035
   for (t in (nyears + 1):(nyears + 2)) {
    Nproc[t] <- log(max(1, lambda*(N[t-1] - H[t-1])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }

})

## Create and compile the NIMBLE model
model <- nimbleModel( code = modelExp,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = "N",
                          control = list(reflective = TRUE),
                          thin = 1,
                         enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
exp_mcmcsamplesPredTests <- runMCMC( mcmc = cMCMC,
                                 nburnin = n.burnin,
                                 niter = n.iter,
                                 nchains = nchains ,
                                 samplesAsCodaMCMC = TRUE,
                                 WAIC = T)
#get sim.lists
exp_mcmcPredTests <- ProcessCodaOutput(exp_mcmcsamplesPredTests$samples,params.omit =c("Ntrueharv","sigmaObs","growth","N"),DIC=F )

#====     3.2.2 exp  Period model ====
  modelPeriod <- nimbleCode({
  # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda ~ dunif(0, 2)
  beta ~ dnorm(0, 1.0E-6)
  N[1] ~ dgamma(1.0E-6, 1.0E-6)
  
  # Process model
  for (t in 2:(nyears)) {
    Nproc[t] <- log(max(1, (lambda + beta*X[t-1])*(N[t-1]- H[t-1])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
  }
  
  # Projected population
  for (t in (nyears + 1):(nyears + 2)) {
    Nproc[t] <- log(max(1, (lambda + beta)*(N[t-1]-H[t-1])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
})


## Create and compile the NIMBLE model
model <- nimbleModel( code = modelPeriod,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                         enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
period_mcmcsamplesPredTests <- runMCMC( mcmc = cMCMC,
                                 nburnin = n.burnin,
                                 niter = n.iter,
                                 nchains = nchains ,
                                 samplesAsCodaMCMC = TRUE,
                                 WAIC = T)
#get sim.lists
period_mcmcPredTests <- ProcessCodaOutput(period_mcmcsamplesPredTests$samples, 
                                          params.omit= c("Ntrueharv","sigmaObs","growth","N"), DIC=F )                                       

#====     3.2.3 Exp  random effect model ====
modelExpRde <- nimbleCode({
 # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda0 ~ dnorm(0, 0.0001) #random growth rate
  lambdaSigma ~ dunif(0,10)
  lambdaTau <- 1/(lambdaSigma^2)
  
  N[1] ~ dgamma(1.0E-6, 1.0E-6)

  # Process model
  for (t in 2:(nyears)) {
    lambdat[t-1] ~ dnorm(lambda0, tau= lambdaTau) # random growth rate
    Nproc[t] <- log(max(1, exp(lambdat[t-1])*(N[t-1]-H[t-1])))#
    N[t] ~ dlnorm(Nproc[t], tauProc)
}

  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
}
  
})
#fit all models
#exp
set.seed(10)
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda0 = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears+2,10,15),
                  N = nimData$Nobs,
                  lambdaSigma= runif(1,0.1,0.2),
                  lambdat = runif(nimConstants$nyears+2,0.1,0.2)
)
nimInits$N <- c(nimInits$N, rep(nimInits$N[nyears-3],NyearsPredic))


NyearsPred <- 2
# on met des NA sur N pour les 2 dernières années. 
thedata1 <- thedata
thedata1$N[(length(thedata1$N)-NyearsPred+1): length(thedata1$N)] <- NA
nimData$Nobs <- thedata1$N
nimInits$Nobs  <- thedata$N
nimInits$Nobs[!is.na(nimData$Nobs)] <- NA
nimConstants$nyears <- nyears
nimParams <- c("N")
nimInits$N <-  thedata$N


## Create and compile the NIMBLE model
model <- nimbleModel( code = modelExpRde,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                          enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
expRde_mcmcsamplesPredTests <- runMCMC( mcmc = cMCMC,
                                        nburnin = n.burnin,
                                        niter = n.iter,
                                        nchains = nchains ,
                                        samplesAsCodaMCMC = TRUE,
                                        WAIC = T)
#get sim.lists
expRde_mcmcPredTests <- ProcessCodaOutput(expRde_mcmcsamplesPredTests$samples, 
                                          params.omit= c("Ntrueharv","sigmaObs","growth","N"), DIC=F )   

#====     3.2.5 plot ====
proj_exp <- exp_mcmcPredTests$sims.list$N %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(medianN = median(value),
            lci = quantile(value, probs = 2.5/100),
            uci = quantile(value, probs = 97.5/100)) %>%
  mutate(an = parse_number(parameter)) %>%
  arrange(an) %>%
  ggplot() + 
  geom_ribbon(aes(x = an, y = medianN, ymin = lci, ymax = uci), fill = "red", alpha = 0.3) + 
  geom_line(aes(x = an, y = medianN), lty = "dashed", color = "red") + 
  #  geom_point(aes(x = an, y = medianN), color = "red") +
  geom_point(data = thedata %>% as_tibble, aes(x = 1:unique(nyears), y = N)) + 
  coord_cartesian(xlim = c(1, nyears), ylim = c(0, 2500)) +
  labs(y = "Effectifs",
       x = "Années",
       title = "Effectifs projetés selon modèle exponentiel",
       subtitle = "avec effectifs observés (points noirs)") + 
  scale_x_continuous(breaks = seq(1, 40, by = 4),
                     labels = seq(1996, 2035, by = 4))



proj_period <- period_mcmcPredTests$sims.list$N %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(medianN = median(value),
            lci = quantile(value, probs = 2.5/100),
            uci = quantile(value, probs = 97.5/100)) %>%
  mutate(an = parse_number(parameter)) %>%
  arrange(an) %>%
  ggplot() + 
  geom_ribbon(aes(x = an, y = medianN, ymin = lci, ymax = uci), fill = "red", alpha = 0.3) + 
  geom_line(aes(x = an, y = medianN), lty = "dashed", color = "red") + 
  geom_point(data = thedata %>% as_tibble, aes(x = 1:unique(nyears), y = N)) + 
  coord_cartesian(xlim = c(1, nyears), ylim = c(0, 2500)) +
  labs(y = "Effectifs",
       x = "Années",
       title = "Effectifs projetés selon modèle exponentiel avec effet aléatoire",
       subtitle = "avec effectifs observés (points noirs)") + 
  scale_x_continuous(breaks = seq(1, 40, by = 4),
                     labels = seq(1996, 2035, by = 4))


proj_expRde <- expRde_mcmcPredTests$sims.list$N %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(medianN = median(value),
            lci = quantile(value, probs = 2.5/100),
            uci = quantile(value, probs = 97.5/100)) %>%
  mutate(an = parse_number(parameter)) %>%
  arrange(an) %>%
  ggplot() + 
  geom_ribbon(aes(x = an, y = medianN, ymin = lci, ymax = uci), fill = "red", alpha = 0.3) + 
  geom_line(aes(x = an, y = medianN), lty = "dashed", color = "red") + 
  geom_point(data = thedata %>% as_tibble, aes(x = 1:unique(nyears), y = N)) + 
  coord_cartesian(xlim = c(1, nyears), ylim = c(0, 2500)) +
  labs(y = "Effectifs",
       x = "Années",
       title = "Effectifs projetés selon modèle frein",
       subtitle = "avec effectifs observés (points noirs)") + 
  scale_x_continuous(breaks = seq(1, 40, by = 4),
                     labels = seq(1996, 2035, by = 4))
library(patchwork)
(proj_exp | proj_period | proj_expRde)  + 
  plot_annotation(title = "Prediction évaluation du modèle")


#====     3.2.6 Calculate relative error ====
relativeErrorExp <- list()
lastYears <- c((nyears-1):nyears)
for(t in 1:length(lastYears) ){
  relativeErrorExp[[t]] <- (exp_mcmcPredTests$sims.list$N[,lastYears[t]]-thedata$N[lastYears[t]])/(thedata$N[lastYears[t]])
}

meanErrorExp <- lapply(relativeErrorExp, mean)
CIErrorExp <- lapply(relativeErrorExp, function(x) quantile(x, probs=c(2.5/100,97.5/100)))
#period
relativeErrorPeriod <- list()
lastYears <- c((nyears-1):nyears)
for(t in 1:length(lastYears) ){
  relativeErrorPeriod[[t]] <- (period_mcmcPredTests$sims.list$N[,lastYears[t]]-thedata$N[lastYears[t]])/(thedata$N[lastYears[t]])
}

meanErrorPeriod <- lapply(relativeErrorPeriod, mean)
CIErrorPeriod <- lapply(relativeErrorPeriod, function(x) quantile(x, probs=c(2.5/100,97.5/100)))

#rde
relativeErrorExpRde <- list()
lastYears <- c((nyears-1):nyears)
for(t in 1:length(lastYears) ){
  relativeErrorExpRde[[t]] <- (expRde_mcmcPredTests$sims.list$N[,lastYears[t]]-thedata$N[lastYears[t]])/(thedata$N[lastYears[t]])
}

meanErrorExpRde <- lapply(relativeErrorExpRde, mean)
CIErrorExpRde <- lapply(relativeErrorExpRde, function(x) quantile(x, probs=c(2.5/100,97.5/100)))


#make a table
tab <- matrix(NA,ncol=3, nrow=3)
row.names(tab) <- c("Constant","Période","Effet Aléatoire") 
colnames(tab) <- c("WAIC","22/23","23/24") 

tab[1,2:3] <- unlist(lapply(1:2, function(x){paste( round(meanErrorExp[[x]], digits=2),"(", round(CIErrorExp[[x]][[1]], digits=2),
                                                 ";", round(CIErrorExp[[x]][[2]], digits=2), ")",sep="")  }))
tab[2,2:3] <- unlist(lapply(1:2, function(x){paste( format(round(meanErrorPeriod[[x]], digits=2),nsmall = 2),"(", round(CIErrorExpRde[[x]][[1]], digits=2),
                                                 ";", round(CIErrorExpRde[[x]][[2]], digits=2), ")",sep="")  }))
tab[3,2:3] <- unlist(lapply(1:2, function(x){paste( round(meanErrorExpRde[[x]], digits=2),"(", round(CIErrorPeriod[[x]][[1]], digits=2),
                                                 ";", round(CIErrorPeriod[[x]][[2]], digits=2), ")",sep="")  }))
#WAIC
tab[1,1] <- round(exp_mcmcsamples$WAIC$WAIC,digits=2)#
tab[2,1] <- round(period_mcmcsamples$WAIC$WAIC,digits=2)#
tab[3,1] <- round(expRde_mcmcsamples$WAIC$WAIC,digits=2)#

write.csv(tab, file = paste(wdFigure,"/WAIC_error.csv",sep=""))


#==== 4. Other plots ====
#====   4.1 raw growth ====
growth <- numeric()
for (i in 2:nyears) {
  growth[i-1] <- thedata$N[i]/thedata$N[i-1]
}
Winter <- paste(substr(year1+1, 3,4),"-",substr(year+1, 3,4),sep="")

###
pdf(file =  paste(wdFigure,"/CroissanceObservée2Periodes.pdf",sep=""),
    width = 10,height = 6)
plot(1996:(1995+length(growth)), growth, type="o", main="",
     xlab="Années", ylab=expression(paste(lambda, " instantané",sep=" ")),ylim=c(0.4,2.5),xaxt="n")
abline(h = 1, col = "gray", lty = 1) # vertical dashed blue line at 2018
axis(1,at=c(1996:(1995+length(growth))),
     labels= Winter[1:length(growth)],las=2,cex.axis=0.8)
abline(v = 2018+0.5, col = "gray", lty = 2) # vertical dashed blue line at 2018

dev.off()


#====   4.2 Prob growth ====
#reduce size of the object by selecting fewer posteriors randomly
set.seed(5)
iterconstant <- sample(dim(exp_mcmc$sims.list$growth)[1],100000) 
growth_mcmc <- exp_mcmc$sims.list$growth[iterconstant, ]
set.seed(5)
iterperiod <- sample(dim(period_mcmc$sims.list$growth)[1],100000) 
growth_mcmPeriod <- period_mcmc$sims.list$growth[iterperiod, ]

Ntrueharv_mcmc <- exp_mcmc$sims.list$Ntrueharv[iterconstant,]
Ntrueharv_mcmcPeriod <- period_mcmc$sims.list$Ntrueharv[iterperiod,]

# Harvest choice - LAMBDA
lower <- 1
upper <- 1.10

p_between <- p_betweenPeriod <- numeric(length(harvestPred))
p_over <- p_overPeriod <- numeric(length(harvestPred))
p_under <- p_underPeriod <- numeric(length(harvestPred))
mean_N <- mean_NPeriod <- numeric(length(harvestPred))
CIL <- CILPeriod <- numeric(length(harvestPred))
CIH <- CIHPeriod <- numeric(length(harvestPred))

# Probability variable x is under a is given by ecdf(x)(a)
for (i in 1:length(harvestPred)) {
  p <- ecdf(growth_mcmc[,i])(upper)
  p_over[i] <- 1 - p
  p_under[i] <- ecdf(growth_mcmc[,i])(lower)
  #Period
  pPeriod <- ecdf(growth_mcmPeriod[,i])(upper)
  p_overPeriod[i] <- 1 - pPeriod
  p_underPeriod[i] <- ecdf(growth_mcmPeriod[,i])(lower)

  #p_under[i] <- mean(growth_mcmc[,i]<1)
mean_N[i] <- mean(Ntrueharv_mcmc[,i])
CIL[i] <- quantile(Ntrueharv_mcmc[,i],probs=c(0.025))
CIH[i] <- quantile(Ntrueharv_mcmc[,i],probs=c(0.975))
p_between[i] <- p - p_under[i]
#Period  
mean_NPeriod[i] <- mean(Ntrueharv_mcmcPeriod[,i])
CILPeriod[i] <- quantile(Ntrueharv_mcmcPeriod[,i],probs=c(0.025))
CIHPeriod[i] <- quantile(Ntrueharv_mcmcPeriod[,i],probs=c(0.975))
p_betweenPeriod[i] <- pPeriod - p_underPeriod[i]

}

#turn harvest to %
HarvestPerc <- array(0,c(dim(exp_mcmc$sims.list$N)[1], length(harvestPred)))
for (i in 1:length(harvestPred)) {
  HarvestPerc[,i] <- harvestPred[i]/exp_mcmc$sims.list$N[,nyears]
}

ptable <- as.data.frame(cbind(harvestPred, colMeans(HarvestPerc), p_under, p_between, p_over, rep(lower, length(harvestPred)), rep(upper, length(harvestPred))))
names(ptable) <- c("Harvest","Harvest %", "P(lambda < 1)", "P(1 < lambda < 1.1)", "P(lambda > 1.1)", "Lower", "Upper")

#====   4.3 Plot growth shooting exp model ====
#exp model 
pdf(file =  paste(wdFigure,"/ProbsversusTirs.pdf",sep=""),
    width = 11,height = 8)
plot(ptable[,1], ptable[,3], ylim=c(0,1), col="red", type="l", 
     xlab="Tirs dérogatoires", ylab="Probabilité", main="")
lines(ptable[,1], ptable[,4], col="grey")
lines(ptable[,1], ptable[,5], col="blue")
legend(15, 0.98, c(expression(paste("P(",lambda, "< 1)",sep="")),
                   expression(paste("P(1 <", lambda, "< 1.1)",sep="")),
                   expression(paste("P(", lambda, "> 1.1)",sep=""))),
                   title= expression(paste( lambda, " entre 2023/24 et 2024/25",sep="")), 
       lty= c(1,1,1), col=c("red", "grey", "blue"), bty="n",cex=1.5)
dev.off()

#
pdf(file =  paste(wdFigure,"/ProbsversusTirsNormalPerc.pdf",sep=""),
    width = 11,height = 8)
plot(ptable[,2]*100, ptable[,3], ylim=c(0,1), col="red", type="l", xlab="Tirs dérogatoires (% de N)", ylab="Probabilité", main="")
lines(ptable[,2]*100, ptable[,4], col="grey")
lines(ptable[,2]*100, ptable[,5], col="blue")
legend(0.05, 0.98, c(expression(paste("P(",lambda, "< 1)",sep="")),
                   expression(paste("P(1 <", lambda, "< 1.1)",sep="")),
                   expression(paste("P(", lambda, "> 1.1)",sep=""))),
       title= expression(paste( lambda, " entre 2023/24 et 2024/25",sep="")), 
       lty= c(1,1,1), col=c("red", "grey", "blue"), bty="n",cex=1.5)
dev.off()



#====   4.3 Table legal mortality vs growth ====
ptable <- as.data.frame(cbind(harvestPred, round(colMeans(HarvestPerc)*100,digits=0), 
                              round(p_under,digits=2), round(p_between,digits=2),round(p_over,digits=2),
                              paste(round(CIL,digits=0),"-",round(CIH,digits=0),sep="")))[seq(6,501,by=10),]
names(ptable) <- c("Harvest", "Harvest %", "P(lambda < 1)",
                   "P(1 < lambda < 1.1)","P(lambda > 1.1)",
                   "Prediction effectifs")
write.csv(ptable[as.numeric(ptable[,1])>50 & as.numeric(ptable[,1])< 400, ],
          file =paste(wdFigure,"/TableTirs.csv",sep="") )





#==== 5. Extra step not in the report; take timing into account in predictions ====
# Number of wolves that can be shot is determined for the year 2025 based on the latest pop size estimates: winter 2023/24.
# However, ideally number of wolves that can shot in 2025 should be defined based on pop size estimates in winter 2024/25.
# this would mean predicting wolf population size in 2024/25 using the number of wolves shot in 2024. 
# This is what is done is Sweden :
# Andrèn H. et al. 2025. Berakningar av beskattning av den skandinaviska vargpopulationen 2026. Rapport till Naturvardsverket och Miljodirektoratet, Norge fran SKANDULV. 41 sidor., 
# this script integrates the known mortality of wolves to predict pop size until 2035. 

#====   5.1 Set data and model  ====
#We add number of wolves that can be shot (192) as data in 2025. 
nimData <- list(Nobs = thedata$N,
                H=thedata$H,
                X = c(rep(0, nyears-8), rep(1, 7), NA),
                harvest= harvestPred
)

nimData$H <- c(nimData$H, 192)

nimConstants <- list(nyears = nyears)
###predict until 2035
NyearsPredic <- 2036 - thedata$year[nyears]
NyearsPred <- nyears + NyearsPredic
nimConstants$NyearsPred <- NyearsPred

set.seed(10)
nimInits <- list( errorObs = runif(1,0,0.2),
                  sigmaProc = runif(1,0.2,1),
                  lambda = runif(1,0.8,1.2),
                  errorObs = runif(1,0,0.2),
                  lambdaNobs = runif(nimConstants$nyears,10,15),
                  beta =  runif(1,-0.1,0.1),
                  N = nimData$Nobs
)

nimInits$N <- c(nimInits$N, rep(nimInits$N[nyears],NyearsPredic))
nimInits$lambdaNobs <- c(nimInits$lambdaNobs, rep(nimInits$lambdaNobs[nyears], NyearsPredic))
nimInits$Ntrueharv <- c(rep(nimInits$lambdaNobs[nyears], length(nimData$harvest)))
nimInits$Nkilled <- c(nimInits$lambdaNobs,rep(nimInits$lambdaNobs[nyears], NyearsPredic))

##model code 
modelExpTiming <- nimbleCode({
  # Priors
  errorObs ~ dunif(0, 0.80)
  sigmaProc ~ dunif(0, 10)
  tauProc <- 1/sigmaProc^2
  lambda ~ dunif(0, 2)
  N[1] ~ dgamma(1.0E-6, 1.0E-6)
  
  # Process model
  for (t in 2:(nyears)) {
    Nproc[t] <- log(max(1, lambda*(N[t-1]- H[t-1])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
  # Observation model
  for (t in 1:nyears) {
    sigmaObs[t] <- errorObs*N[t]
    shapeObs[t] <- N[t]*N[t]/(sigmaObs[t]*sigmaObs[t])
    rateObs[t] <- N[t]/(sigmaObs[t]*sigmaObs[t])
    lambdaNobs[t] ~ dgamma(shapeObs[t], rateObs[t])
    Nobs[t] ~ dpois(lambdaNobs[t])
  }
  
  # Projected population until 2035
  for (t in (nyears + 3):NyearsPred) {
    #here we fix the rate to 19%
    Nkilled[t] <- N[t-1] *0.19
    Nproc[t] <- log(max(1, lambda*(N[t-1]-Nkilled[t])))
    N[t] ~ dlnorm(Nproc[t], tauProc)
  }
  
  # at Nyears+1, the number of wolves dead is known. We can use it 
  Nproc[nyears + 1] <- log(max(1, lambda*(N[nyears]-H[nyears])))
  N[nyears + 1] ~ dlnorm(Nproc[nyears + 1], tauProc)
  # at Nyears+2, the number of wolves dead is also known. We can use it 
  Nproc[nyears + 2] <- log(max(1, lambda*(N[nyears+1]-H[nyears+1])))
  N[nyears + 2] ~ dlnorm(Nproc[nyears + 2], tauProc)
  
  ##################################
  # Harvested population
  ##################################
  #Instead here we make the prediction of harvest on the predicted pop size at nyears+2 given different rate of harvest.
  #(even if we know the number of wolves that will likely be shot (192))
  for (i in 1:length(harvest)) {
    Nprocharv[i] <- log(max(1, lambda*(Nproc[nyears + 2] - harvest[i]))) # Because of index
    Ntrueharv[i] ~ dlnorm(Nprocharv[i], tauProc)
    growth[i] <- Ntrueharv[i]/ N[nyears + 2]
  }
})

###
nimParams <- c("lambda", "sigmaProc", "sigmaObs", "N", "tauProc","Ntrueharv", "growth")

## Create and compile the NIMBLE model
model <- nimbleModel( code = modelExpTiming,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)
## Check the initial log-likelihood 
model$calculate()
#compile and build model
cmodel <- compileNimble(model)
MCMCconf <- configureMCMC(model = model,
                          monitors  = nimParams,
                          control = list(reflective = TRUE),
                          thin = 1,
                          enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
## Run MCMC
set.seed(5)
exp_mcmcsamplesTiming <- runMCMC( mcmc = cMCMC,
                                  nburnin = n.burnin,
                                  niter = n.iter,
                                  nchains = nchains ,
                                  samplesAsCodaMCMC = TRUE,
                                  WAIC = T)
#get sim.lists
exp_mcmcTiming <- ProcessCodaOutput(exp_mcmcsamplesTiming$samples,
                                    params.omit= c("Ntrueharv","sigmaObs","growth","N"), DIC=F )                                       

#====   4.2 Plots  ====
  dat <- exp_mcmcTiming$sims.list$N %>%
  as_tibble() %>%
  pivot_longer(cols = everything(),  values_to = "value", names_to = "parameter") %>%
  group_by(parameter) %>%
  summarize(medianN = median(value),
            lci = quantile(value, probs = 2.5/100),
            uci = quantile(value, probs = 97.5/100),
            lci50 = quantile(value, probs = 25/100),
            uci50 = quantile(value, probs = 75/100),
            lci75 = quantile(value, probs = 12.5/100),
            uci75 = quantile(value, probs = 87.5/100),
  ) %>%
  mutate(an = parse_number(parameter)) %>%
  arrange(an) 


#
proj_expTiming <- 
  ggplot() + 
  geom_ribbon(data= dat, aes(x = an, y = medianN, ymin = lci, ymax = uci), 
              fill = grey(0.6), alpha = 0.4) + 
  geom_ribbon(data= dat, aes(x = an, y = medianN, ymin = lci75, ymax = uci75), 
              fill = grey(0.4), alpha = 0.4) +
  geom_ribbon(data= dat, aes(x = an, y = medianN, ymin = lci50, ymax = uci50), 
              fill = grey(0.1), alpha = 0.4) + 
  
  # geom_ribbon(data= dat[nyears: nrow(dat),], aes(x = an, y = medianN, ymin = lci50, ymax = uci50), 
  #             fill = grey(0.3), alpha = 0.6) + 
  #geom_line(data= dat[1:nyears+0,], aes(x = an, y = medianN), lty = "dashed", color = "red") + 
  #  geom_point(aes(x = an, y = medianN), color = "red") +
  geom_point(data = thedata %>% as_tibble, aes(x = 1:unique(nyears), y = N)) + 
  coord_cartesian(xlim = c(1, 40), ylim = c(0, 4200)) +
  labs(y = "Effectifs",
       x = "Hivers")+
  #title = "Effectifs projetes selon modèle exponentiel",
  #subtitle = "avec effectifs observes (points noirs)") + 
  scale_x_continuous(breaks = seq(1, 41, by = 4),
                     labels = Winter[seq(1, 41, by = 4)])

#Incorporating known mortality does not add much value to the predictions at 2035
(proj_exp2035Line5075 | proj_expTiming )  + 
  plot_annotation(title = "Prediction évaluation du modèle")




