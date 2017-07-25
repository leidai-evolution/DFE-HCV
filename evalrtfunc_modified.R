# These are the supporting functions for evalrt.R and evalrtmc.R
# which implement methods described in Beisel et. al.
# Testing the extreme value domain of attraction for distributions of beneficial fitness effects." Genetics (2007)

#updated: Lei Dai, 07/11/2017

dgpd<-function (x, tau = 1, kappa = 0) 
{
  if (tau <= 0) stop("tau must be positive")
  if (kappa < 0 && x <= 0 && x > -(tau/kappa)) d<- -Inf
  else {
    if(kappa == 0) d<- exp(-x/tau)/tau
    else d<- ((1 + kappa*x/tau)^(-(kappa+1)/kappa))/tau
  }
  d
}

gpdnegloglike<-function (theta,y) #y should be a parameter here
{
  tau<-exp(theta[1])
  kappa<-exp(theta[2])-1
  n<-length(y)
  if (kappa == 0) {
    l <- -(n)*log(tau)-(1/tau)*sum(y)
  }
  else {
    if(min(1+(kappa/tau)*y)<=0) l<- -Inf
    else{l <- -(n)*log(tau)-((kappa+1)/kappa)*sum(log(1+(kappa/tau)*y))}	
  } 
  -l
}

expnegloglike<-function (tau,y) #y should be a parameter here
{
  tau<-exp(tau)
  n<-length(y)
  l <- -(n)*log(tau)-(1/tau)*sum(y) 
  -l
}

rgpd <- function (n, tau = 1, kappa = 0) {
  if (tau <= 0) stop("tau must be positive")
  if (abs(kappa) < .001) return(tau * rexp(n))
  else return(tau * (runif(n)^(-kappa) - 1)/kappa)
}

gpdneglogmclike<-function(theta, y, x, tauhat, kappahat) {  #y,x,tauhat,kappahat should be a parameter here
  tau<-exp(theta[1])
  kappa<-exp(theta[2])-1
  n<-length(y)
  l<-dnorm(y,x,sigma)*dgpd(x,tau,kappa)/dgpd(x,exp(tauhat),exp(kappahat)-1)
  like<-sum(log(apply(l,1,FUN=sum)/n))
  -like
}

expneglogmclike<-function(tau, y, x, tauhat, kappahat) { #y,x,tauhat,kappahat should be a parameter here
  kappa<-0
  tau<-exp(tau)
  n<-length(y)
  l<-dnorm(y,x,sigma)*dgpd(x,tau,kappa)/dgpd(x,exp(tauhat),exp(kappahat)-1)
  like<-sum(log(apply(l,1,FUN=sum)/n))
  -like
}
