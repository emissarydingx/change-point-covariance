library(MASS)
library(compiler)
library(doParallel)
library(foreach)
library(mvtnorm)
library(ContaminatedMixt)
# library(distrEllipse)

# #######################################
# 1. bootstrap_gen
# n: n samples
# d: d dimensions
# nBoot: a number represents the number of sampling
# cp_indx: the change point index (after that the covariance matrices begins changing)
# Sigmas: a list that contains covariance matrix of data
bootstrap_gen<-function(n,d,cp_indx,family,para){
  Sigma1=para$Sigma1
  Sigma2=para$Sigma2
  #normal distribution
  if (family=='normal'){
    X1=mvrnorm(n = cp_indx, mu=rep(0,d), Sigma=Sigma1)
    X2=mvrnorm(n = n-cp_indx, mu=rep(0,d), Sigma=Sigma2)
  }
  #t distribution
  if (family=='t'){
    X1=rmvt(cp_indx, sigma = Sigma1, df = para$df)
    X2=rmvt(n-cp_indx, sigma = Sigma2, df = para$df)
  }
  # mixed distribution
  if (family=='mix'){
    X1=rCN(n = cp_indx, mu = rep(0,d), Sigma1, alpha = 0.8, eta = 1.01)
    X2=rCN(n = n-cp_indx, mu = rep(0,d), Sigma2, alpha = 0.8, eta = 1.01)
  }
  #data matrix
  X=rbind(X1,X2)
  return(X)
}


# #######################################
# 2. bootstrap_fit
# X: data matrix (n by p)
# it: it th observation
# nBoot: a number represents the number of sampling
# b: bandwidth

dev_Sig_at_t_with_e<-function(s,X_star, e){
  n=dim(X_star)[1]
  d=dim(X_star)[2]
  s1=crossprod(X_star[(1:s),]*drop(e[(1:s)]),X_star[(1:s),])-
    sum(e[1:s])/s*crossprod(X_star[(1:s),],X_star[(1:s),])
  s2=crossprod(X_star[((s+1):n),]*drop(e[((s+1):n)]),X_star[((s+1):n),])-
    sum(e[(s+1):n])/(n-s)*crossprod(X_star[((s+1):n),],X_star[((s+1):n),])
  dev=max(abs(sqrt((n-s)/(n*s))*s1-sqrt(s/(n*(n-s)))*s2))
  return(dev)
}

bootstrap_fit<-function(X,nBoot,s0){
  n=dim(X)[1]
  d=dim(X)[2]
  #generate iid e_i: nBoot_e n by nBoot
  # nBoot_e=mvrnorm(n = n, mu=rep(0,nBoot), Sigma=diag(rep(1,nBoot)))
  nBoot_e=sapply(1:nBoot,FUN=function(x){
    rnorm(n, mean = 0, sd = 1)
  })
  #bootstraping
  # X_bootstrap: n*d*nBoot
  X_bootstrap=sapply(1:nBoot,FUN=function(x){
    indx_tmp=sample(1:n,replace=TRUE)
    return(X[indx_tmp,])
  },simplify="array")  
  tmp=sapply(1:nBoot,FUN=function(step){
    sapply(s0:(n-s0),FUN=function(s){
      X_star=X_bootstrap[,,step]
      e=nBoot_e[,step]
      return(dev_Sig_at_t_with_e(s,X_star,e))
    })
  },simplify="array")
  
  W_all=as.vector(apply(tmp,MARGIN=2,max))
  return(W_all)
}

dev_Sig_at_t_without_e<-function(s,X){
  n=dim(X)[1]
  d=dim(X)[2]
  s1=crossprod(X[1:s,],X[1:s,])
  s2=crossprod(X[(s+1):n,],X[(s+1):n,])
  dev=max(abs(sqrt((n-s)/(n*s))*s1-sqrt(s/(n*(n-s)))*s2))
  # dev=norm(s1-s2,type='I')/hnb
  return(dev)
}

T_statistic<-function(X,s0){
  n=dim(X)[1]
  d=dim(X)[2]
  tmp=sapply(s0:(n-s0),FUN=function(s){
    return(dev_Sig_at_t_without_e(s,X))
  })
  return(max(tmp))
}








