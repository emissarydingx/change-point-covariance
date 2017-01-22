rm(list=ls())
library(MASS)
library(doParallel)
library(foreach)
library(mvtnorm)
library(ContaminatedMixt)
# library(distrEllipse)

# setwd("C:\\Users\\Xin\\OneDrive\\research\\change point covariance\\code\\simulation")
# setwd('/home/xin/OneDrive/Research/change point covariance/code/simulation')
# setwd("/home/xin/Dropbox/Research/change point covariance/code/simulation")
setwd("C:\\Users\\Xin\\Dropbox\\Research\\change point covariance\\code\\simulation")
source("change_point_bootstrap_lib_old.R")
# -------------------------------------------
# Simulation setup
# -------------------------------------------
ncores=8
n=100
d=40
nBoot=200
nSim=100

#bandwidth fitting
# b=0.8*n^(-1/5)
b=0.3

Sigma1=diag(rep(10,d))
Sigma2=diag(rep(10,d))


cp_indx=50#the index of the change point


para=list(Sigma1=Sigma1,Sigma2=Sigma2,
          df=10)

alpha_hat_arr=array(0,c(nSim,25))
alpha=seq(0.01,0.99,length.out=25)#range of probability

ptm <- proc.time()
cl <- makeCluster(ncores)
registerDoParallel(cl)
for (nn in 1:nSim){
  print(nn)
  set.seed(nn)
  #data generation
  X=bootstrap_gen(n,d,cp_indx,family='t',para=para)
  #bootstrap
  W=bootstrap_fit(X,nBoot,b,ncores)
  W_sort=sort(W)
  #Empirical cdf of W
  ecdf_W=ecdf(W)
  prob_W=ecdf_W(W_sort)#discontinuity points of W's ecdf
  prob_W_interval=c(0,prob_W[-length(prob_W)])
  #sample quantile (prob 0 represents smallest sample while prob 1 represents largest sample)
  sq=sapply(1:length(alpha),FUN=function(x){
    W_sort[sum(alpha[x]>prob_W_interval)]
  })
  # sq=quantile(W,probs=alpha)
  T_n=T_statistic(X,b,ncores)
  alpha_hat_arr[nn,]=sapply(1:length(sq),FUN=function(x){
    T_n<=sq[x]
  })
}
stopCluster(cl)
t=proc.time() - ptm
# alpha_hat=as.vector(apply(alpha_hat_arr,MARGIN=2,mean))
alpha_hat=colMeans(alpha_hat_arr)

plot(alpha,alpha_hat,xlab='alpha',ylab='Bootstrap approximation',type='b',xlim=c(0, 1),ylim=c(0,1))
abline(a=0,b=1)

# save.image("change_point_detection.Rdata")



