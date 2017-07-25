# this script performs:
# 1) Maximum Likelihood Estimate of the full model (Generalized Pareto Distribution) 
# and reduced model (exponential dist)
# 2) Likelihood Ratio Test: return p-value
# updated: 09/21/2015, LD

rm(list=ls())
source("gpdmle.R")

#####MAIN#####
# input: HCV NS5A fitness data, output of analysis_RF2.m
inputfile <- "data/fitness_singleaa.txt"
inputfile2 <- "data/fitness_pointmut.txt"
inputfile_synonymous <- "data/fitness_synonymous.txt"
# non-synonymous mutations
table_singleaa <- read.table(inputfile,header = TRUE,sep = "\t")
table_pointmut <- read.table(inputfile2,header = TRUE,sep = "\t")
# synonymous mutations: used to determine the threshold for beneficial mutations
# silent <- scan(inputfile2,sep = "\t")
table_silent <- read.table(inputfile_synonymous,header = TRUE,sep = "\t")
#fitness_silent <- table_silent[,2]
#s=log(Relative Fitness)
temp <- table_silent[,2] 
fitness_silent <- log(temp[temp>0])

#fit the distribution of beneficial mutations####
#single aa mutations
output <- matrix(,nrow=4,ncol=5)
num_ben <- matrix(,nrow=4,ncol=1)
frac_ben <- matrix(,nrow=4,ncol=1)

for (i in 2:ncol(table_singleaa)){
  #identify benefical mutations for each condition (0,10,40,100pM)
  #fitness <- table_singleaa[,i]
  #log
  temp <- table_singleaa[,i]
  fitness <- log(temp[temp>0])
  #threshold <- 1+2*sd(fitness_silent) #specify fitness threshold for beneficial mutations
  threshold <- 2*sd(fitness_silent) #specify fitness threshold for beneficial mutations
  fitness_ben <- fitness[fitness>threshold]
  num_ben[i-1] <- length(fitness_ben)
  frac_ben[i-1] <- length(fitness_ben)/length(fitness) 
  
  fitness_min <- min(fitness_ben)
  s <- fitness_ben-fitness_min #shift the distribution relative to the smallest observed selection coefficient
  # set number of bootstrap replicates to perform
  bootsize <- 10000
  output[i-1,]  <- gpdmle(s,bootsize)  
}

#output####
#outputfile <- "data/GPDfit_singleaa.txt"
outputfile <- "data/GPDfit_singleaa_log.txt"
header_par <- c('tau_full','kappa_full','tau_reduced','kappa_reduced','pvalue','num_ben','frac_ben')
write.table(cbind(output,num_ben,frac_ben), file = outputfile, sep = "\t", col.names = header_par, row.names =FALSE)

#single nt mutations####
num_ben_pointmut <- matrix(,nrow=4,ncol=1)
frac_ben_pointmut <- matrix(,nrow=4,ncol=1)
output_pointmut <- matrix(,nrow=4,ncol=5)
for (i in 2:ncol(table_pointmut)){
  #identify benefical mutations for each condition (0,10,40,100pM)
  fitness <- table_pointmut[,i]
  threshold <- 1+2*sd(fitness_silent) #specify fitness threshold for beneficial mutations
  fitness_ben <- fitness[fitness>threshold]
  num_ben_pointmut[i-1] <- length(fitness_ben)
  frac_ben_pointmut[i-1] <- length(fitness_ben)/length(fitness)
  
  fitness_min <- min(fitness_ben)
  s <- fitness_ben-fitness_min #shift the distribution relative to the smallest observed selection coefficient
  # set number of bootstrap replicates to perform
  bootsize <- 10000
  output_pointmut[i-1,]  <- gpdmle(s,bootsize)  
}

#output_pointmut####
#fit using Relative Fitness
#outputfile <- "data/GPDfit_pointmut.txt"
#fit using log(Relative Fitness)
outputfile <- "data/GPDfit_pointmut_log.txt"
header_par <- c('tau_full','kappa_full','tau_reduced','kappa_reduced','pvalue','num_ben','frac_ben')
write.table(cbind(output_pointmut,num_ben_pointmut,frac_ben_pointmut), file = outputfile, sep = "\t", col.names = header_par, row.names =FALSE)

