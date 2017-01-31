rm(list=ls())
library(MASS)
library(doParallel)
library(foreach)
library(mvtnorm)
library(ContaminatedMixt)
# library(distrEllipse)

# setwd("/home/xin/Dropbox/Research/change point covariance/code/simulation")
setwd("C:\\Users\\Xin\\Dropbox\\Research\\change point covariance\\code\\simulation")
source("change_point_bootstrap_lib.R")
# -------------------------------------------
# Simulation setup
# -------------------------------------------
ncores=8
n=100
d=20
nBoot=300
nSim=100
s0=5
family='mix'

filename=paste0('cp_n_',n,'_d_',d,'_nBoot_',nBoot,'_nSim_',nSim,'_s0_',s0,'_family_',family,'.Rdata',collapse = NULL)

Sigma1=diag(rep(1,d))
Sigma2=diag(rep(1,d))


cp_indx=50#the index of the change point


para=list(Sigma1=Sigma1,Sigma2=Sigma2,
          df=5)

alpha_hat_arr=array(0,c(nSim,25))
alpha=seq(0.01,0.99,length.out=25)#range of probability

#data generation
X_all=array(0,c(n,d,nSim))
for (nn in 1:nSim){
  set.seed(nn)
  X_all[,,nn]=bootstrap_gen(n,d,cp_indx,family=family,para=para)
}


ptm <- proc.time()
cl <- makeCluster(ncores)
registerDoParallel(cl)
alpha_hat_arr<-foreach(nn=1:nSim,.combine='rbind',
                       .packages=c('MASS','mvtnorm','ContaminatedMixt'))%dopar%{
  X=X_all[,,nn]
  #bootstrap
  W=bootstrap_fit(X,nBoot,s0)
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
  T_n=T_statistic(X,s0)
  array(sapply(1:length(sq),FUN=function(x){
    T_n<=sq[x]
  }),c(1,length(sq)))
}
stopCluster(cl)
t=proc.time() - ptm
alpha_hat=as.vector(apply(alpha_hat_arr,MARGIN=2,mean))
alpha_hat=colMeans(alpha_hat_arr)

plot(alpha,alpha_hat,xlab='alpha',ylab='Bootstrap approximation',type='b',xlim=c(0, 1),ylim=c(0,1))
abline(a=0,b=1)

# save.image("change_point_detection.Rdata")

save.image(filename)



