# updated: 09/21/2015, LD

gpdmle <- function (s,bootsize){
  source("evalrtfunc_modified.R")
  
  y<-sort(s)
  # Compute sample size n
  n<-length(y)
  cat("n=",n,"\n")
  print(y)
  
  #MLE
  reducedtau<-mean(y)
  reduced<-expnegloglike(log(reducedtau),y)
  full<-optim(c(log(reducedtau),0),gpdnegloglike,y=y)  
  lrtts<-(-2*(full$value-reduced))
  cat("Full Model: ")
  cat("tauhat = ",exp(full$par[1]),"\tkappahat = ",exp(full$par[2])-1,"\n")
  cat("Reduced Model: ")
  cat("tauhat = ",reducedtau,"\tkappahat = ",0,"\n")
  cat("l(full)=\t",-full$value,"\n")
  cat("l(reduced)=\t",-reduced,"\n")
  lrtts<-(-2*(full$value-reduced))
  cat("LRT\t\t=\t",lrtts,"\n\n")
  
  #tau, kappa: full model and reduced model (exponential)
  #full: 1st row; reduced: 2nd row
  MLE_full <-c(exp(full$par[1]),exp(full$par[2])-1)
  MLE_reduced <-c(reducedtau,0)
  #  print(MLE_full)
  #  print(MLE_reduced)
  
  #LRT
  #Generate Parametric Bootstrap of Test Statistic
  teststat <- rep(0,bootsize)
  for(i in 1:bootsize) {
    # Generate data set of size n
    y_temp<-rgpd(n,reducedtau,0) #random variates from the GEV distribution
    y_temp<-sort(y_temp)
    # Perform Likelihood Ratio Test
    boottau<-mean(y_temp)
    reduced<-expnegloglike(log(boottau),y_temp)
    full<-optim(c(log(boottau),0),gpdnegloglike,y=y_temp)
    teststat[i]<- (-2*(full$value-reduced))
  }
  #hist(teststat)
  f <- ecdf(teststat)
  pvalue<-1-f(lrtts)
  cat("p-value\t\t=\t",pvalue,"\n")
  
  return(c(MLE_full,MLE_reduced,pvalue))
  
}