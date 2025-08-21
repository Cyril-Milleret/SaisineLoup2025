###############################################
# WOLF MATRIX POPULATION MODELS
###############################################
#SCRIPT ADAPATED FROM Marescot et al 2012 Reducing matrix population models with application to social animal species. Ecological Modelling
#   232 :  91-96. https://doi.org/10.1016/j.ecolmodel.2012.02.017

rm(list=ls())

###############################################
# Demographic parameters
###############################################
MCiter<-1000000
MCiter<-50000
# Reproduction
litter_size <- runif(MCiter,4,9)
f <- litter_size/2  
# Dispersal and installation
p_es <- runif(MCiter, 0.3, 0.7)
p_di <- runif(MCiter, 0.1, 0.4)

# Mortality
mortality <- seq(0.05, 0.8, 0.025)
phi<-matrix(0,nrow=MCiter, ncol=length(mortality))
for (d in 1:length(mortality))
{
  phi[,d] <- 1 - mortality[d]
  
  phi[,d] <- rnorm(MCiter,phi[,d],0.1)
  index1 <- which(phi[,d]<0)
  index2 <- which(phi[,d]>1)
  
  while (length(index1) > 0) 
  {    phi[index1,d] <- rnorm(length(index1), 1 - mortality[d], 0.1)
  index1 <- which(phi[,d]<0)
  }    
  while(length(index2) >0 )
  {
    phi[index2,d] <- rnorm(length(index2), 1 - mortality[d], 0.1) 
    index2 <- which(phi[,d]>1)
  }
  
}
###############################################
# MODEL WITH FOUR CLASSES: J, S, D, A
###############################################

###############################################
# Loop on mortality with matrix model
###############################################

NStages <- 4

pgr4 <- rep(NA, length(mortality))
lower <- rep(NA, length(mortality))
upper <- rep(NA, length(mortality))
lambda<-matrix(0,nrow=MCiter, ncol=length(mortality))

for(d in 1:length(mortality)) {
  popvector <- rep(0,times=NStages)
  
  M <- array(data=0, dim=c(NStages, NStages, MCiter)) 
  
  
  phi_j <- phi[,d]
  phi_s <- phi[,d]
  phi_a <- phi[,d]
  phi_d <- phi[,d]
  
  
  M[1,4,]= f*phi_a
  M[2,1,]= phi_j*p_di
  M[2,3,]= phi_s
  M[3,1,]= phi_j*(1-p_di)
  M[4,2,]= phi_d*p_es
  M[4,4,]= phi_a
  
  for(i in 1:MCiter)
  {
    Mat=matrix(nrow=NStages, ncol=NStages, data=M[,,i]) 
    rightV <-eigen(Mat)                     
    lambda[i,d] <- log(Re(rightV$values[1]))
  }
  pgr4[d] <- median(lambda[,d])
  lower[d] <- quantile(lambda[,d],probs=2.5/100)
  upper[d] <- quantile(lambda[,d], probs=97.5/100)
}

###############################################
# Data (all) and graph
###############################################

fullerdata <- c(0.68,0.46,0.45,0.31,0.36,0.33,0.36,0.28,0.15)
fullerdata <- rev(fullerdata) 
marescotdata <- 0.161
maruccodata <- 0.18
mech_boitani <- c(0.68, 0.58, 0.56, 0.45, 0.45, 0.42, 0.31, 0.60, 0.34, 0.21, 0.37, 0.36, 0.33, 0.27, 0.35, 0.28, 0.18, 0.15, 0.16)
#
mort <- c(maruccodata, marescotdata, mech_boitani, fullerdata)

fullerresult <- c(-0.92,-0.37,-0.13,-0.08,0.02,0.06,0.1,0.17,0.19)
fullerresult <- rev(fullerresult)
marescotresult <- 0.271
maruccoresult <- (32/12)^(1/7)-1 
mech_boitaniresult <- c(-0.92, -0.27,-0.15,-0.13,-0.13, -0.12,-0.08,-0.03,-0.05,0.01,0.01,0.02,0.03,0.18,0.06,0.12,0.15,0.19,0.4)

result <- c(maruccoresult, marescotresult, mech_boitaniresult, fullerresult)


#plot with growth as an exponential scale
pgr4exp <- exp(pgr4)
resultexp <- exp(result)
upperexp <- exp(upper)
lowerexp <- exp(lower)

plot(pgr4exp ~ mortality, main="", ylab='Taux de croissance',
     xlab='Taux de mortalité', ylim=c(0,2), col='red', type="o")
##add extra points. Because survival was age-sepcific and we did not know the % of individuals in age class, 
# we assumed that overall mortality was equal to the mean young and adult individusl
#-	SCANDINAVIA (Mileret 2014)  L=0.91 Smark=0.61,Sother= 0.46 > Smean= 0,535
#-	ITALY (Marucco et al 2009) L=1.04 Sad=0,82 Sj=0.24 > Smean=0,53   
#-	GERMANY (planilo et al  2024) : L=1,31 Sad 0,88 Ssubj=0,74 > Smean=0,81
mort1  <- c(mort, 1-0.535, 1-0.53, 1-0.81 )
resultexp1 <- c(resultexp, 0.91, 1.04, 1.31 )
points(mort1, resultexp1, pch=20)
abline(h=1, lty=2,col=grey(0.7))
#plot France 
points(0.38,1.11, pch=20)#, col="blue")
growthIC <- c(1.08,1.22)
mortIC <- 1-c(0.56,0.69)
lines(mortality,  upperexp, lty="dashed", col="red")
lines(mortality,  lowerexp, lty="dashed", col="red")

#flags were added after printing the figure 
