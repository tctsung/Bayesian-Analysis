library(MCMCpack)
library(coda)
#### Social mobility data
load("socmob.RData") 

yincc<-match(socmob$INC,sort(unique(socmob$INC)))
ydegr<-socmob$DEGREE+1
ychild<-socmob$CHILD
ypdeg<-1*(socmob$PDEG>2)

#####
X<-cbind(ychild,ypdeg,ychild*ypdeg)
y<-ydegr
data<-na.omit(cbind(y,1,X)) # remove missing data
y<-data[,1]
X<-data[,-1]

ranks<-match(y,sort(unique(y))) 
uranks<-sort(unique(ranks))
K<-length(uranks)

n<-dim(X)[1] ; p<-dim(X)[2]
XX<-t(X)%*%X  

#### Ordinal probit regression #############
## setup
beta.0<-rep(0,p) # prior mean on beta
Sigma.0<-n*solve(XX) # prior covariance matrix on beta
iSigma.0<-solve(Sigma.0)
mu.0<-rep(0,K-1) # prior mean on g
gamma.0<-rep(10,K-1) #prior SD on g

beta<-rep(0,p) #initial values

Data<-data.frame(ydegr=y,ychild=X[,2],ypdeg=X[,3],ychild.ypdeg=X[,4])

## MCMC
S<-26000
# An intercept is needed and assumed in MCMCoprobit()
model<-MCMCoprobit(ydegr~ychild+ypdeg+ychild.ypdeg, data=Data,
                  burnin=0, 
                  mcmc=S, 
                  seed=0, # random seed for random number generators
                  tune=0.2/K,# The tuning parameter for the Metropolis-Hastings step. 
                  # Default of NA corresponds to a choice of 
                  # 0.05 divided by the number of categories in the response variable.
                  beta.start=beta, # initial values of beta
                  b0=beta.0, # prior mean of beta
                  B0=iSigma.0, # prior precision matrix (inverse covariance matrix) of beta
                  a0=mu.0, # prior mean of gamma. Here, gamma-intercept=g. Note that gamma1=0 is fixed.
                  A0=diag(1/gamma.0^2), # prior precision matrix of gamma
                  verbose=1)
# The following warning doesn't matter
# Warning message:
#  In if (is.na(beta.start)) { :
#      the condition has length > 1 and only the first element will be used


attributes(model) # burnin samples are discarded.
summary(model)
colnames(model)
p=dim(model)[2]
mcmc.pars<-matrix(0,S,p)
colnames(mcmc.pars)<-colnames(model)
for(i in 1:p){
  mcmc.pars[,i]<- as.numeric(model[,i]) # change to numeric, otherwise the index uses the iteration # even after thinning
}



# Traceplot
par(mfrow=c(p,1))
for(j in 1:p){# plot the first 2000 iterations
  plot(mcmc.pars[1:2000,j],ylab=colnames(mcmc.pars)[j],xlab="Iteration",main=paste("Trace plot for ",colnames(mcmc.pars)[j],sep=""),type="l",lwd=1,lty=1)
}


#ACF plot
burnin<-1000
par(mfrow=c(p-K+1,1))
for(j in 1:(p-K+1)){
  acf(mcmc.pars[(burnin+1):S,j],main=paste("ACF for ",colnames(mcmc.pars)[j],sep=""),lag.max=100)
}


thin<-25
par(mfrow=c(p-K+1,1))
for(j in 1:(p-K+1)){# after thinning
  acf(mcmc.pars[seq(burnin+1,S,by=thin),j],main=paste("ACF for ",colnames(mcmc.pars)[j],sep=""),lag.max=100)
}

length(mcmc.pars[seq(burnin+1,S,by=thin),1])
for(j in 1:(p-K+1)){
print(effectiveSize(mcmc.pars[seq(burnin+1,S,by=thin),j]))
}


mcmc.pars.pMean<-apply(mcmc.pars[seq(burnin+1,S,by=thin),],2,mean) # posterior mean of parameters
mcmc.pars.pMean
mcmc.pars.pSD<-apply(mcmc.pars[seq(burnin+1,S,by=thin),],2,sd) # posterior mean of parameters
mcmc.pars.pSD
# Note that gamma1=0 is fixed by MCMCpack, and gamma-intercept=g
# Compare with Lecture11_code_1.R

########### for binary case ##############
y2=as.numeric(Data$ydegr>3)
Data2<-data.frame(ydegr=y2,Data[,2:4])

###* probit regression ####
# MCMCprobit uses the Gibbs sampler
model2<-MCMCprobit(ydegr~ychild+ypdeg+ychild.ypdeg, data=Data2,
                   burnin=0, 
                   mcmc=S, 
                   seed=0, # random seed for random number generators
                   beta.start=beta, # initial values of beta
                   b0=beta.0, # prior mean of beta
                   B0=(iSigma.0+t(iSigma.0))/2 # prior precision matrix (inverse covariance matrix) of beta
                  # B0 needs a symmetric matrix. 
                  # But iSigma.0 is not exactly symmetric due to computational precision.
                  # So, we use (iSigma.0+t(iSigma.0))/2 which is symmetric.
                  )


attributes(model2) # burnin samples are discarded.
summary(model2)
colnames(model2)
p=dim(model2)[2]
mcmc2.pars<-matrix(0,S,p)
colnames(mcmc2.pars)<-colnames(model2)
for(i in 1:p){
  mcmc2.pars[,i]<- as.numeric(model2[,i]) # change to numeric, otherwise the index uses the iteration # even after thinning
}



# Traceplot
par(mfrow=c(p,1))
for(j in 1:p){# plot the first 2000 iterations
  plot(mcmc2.pars[1:2000,j],ylab=colnames(mcmc2.pars)[j],xlab="Iteration",main=paste("Trace plot for ",colnames(mcmc2.pars)[j],sep=""),type="l",lwd=1,lty=1)
}


#ACF plot
burnin<-1000
par(mfrow=c(p,1))
for(j in 1:p){
  acf(mcmc2.pars[(burnin+1):S,j],main=paste("ACF for ",colnames(mcmc2.pars)[j],sep=""),lag.max=100)
}


thin<-5
par(mfrow=c(p,1))
for(j in 1:p){# after thinning
  acf(mcmc2.pars[seq(burnin+1,S,by=thin),j],main=paste("ACF for ",colnames(mcmc2.pars)[j],sep=""),lag.max=100)
}

length(mcmc2.pars[seq(burnin+1,S,by=thin),1])
for(j in 1:p){
  print(effectiveSize(mcmc2.pars[seq(burnin+1,S,by=thin),j]))
}


mcmc2.pars.pMean<-apply(mcmc2.pars[seq(burnin+1,S,by=thin),],2,mean) # posterior mean of parameters
mcmc2.pars.pMean
mcmc2.pars.pSD<-apply(mcmc2.pars[seq(burnin+1,S,by=thin),],2,sd) # posterior mean of parameters
mcmc2.pars.pSD


###* logistic regression ####
# MCMClogit use a random walk Metropolis algorithm
model3<-MCMClogit(ydegr~ychild+ypdeg+ychild.ypdeg, data=Data2,
                  burnin=0, 
                  mcmc=S, #
                  seed=0, # random seed for random number generators
                  tune= 1.1, # tuning parameter for the Metropolis algorithm. 
                  # The default value is 1.1
                  beta.start=beta, # initial values of beta
                  b0=beta.0, # prior mean of beta
                  B0=(iSigma.0+t(iSigma.0))/2, # prior precision matrix (inverse covariance matrix) of beta
                  # B0 needs a symmetric matrix. 
                  # But iSigma.0 is not exactly symmetric due to computational precision.
                  # So, we use (iSigma.0+t(iSigma.0))/2 which is symmetric.
                  verbose=1
)

attributes(model3) # burnin samples are discarded.
summary(model3)
colnames(model3)
p=dim(model3)[2]
mcmc3.pars<-matrix(0,S,p)
colnames(mcmc3.pars)<-colnames(model3)
for(i in 1:p){
  mcmc3.pars[,i]<- as.numeric(model3[,i]) # change to numeric, otherwise the index uses the iteration # even after thinning
}



# Traceplot
par(mfrow=c(p,1))
for(j in 1:p){# plot the first 2000 iterations
  plot(mcmc3.pars[1:2000,j],ylab=colnames(mcmc3.pars)[j],xlab="Iteration",main=paste("Trace plot for ",colnames(mcmc3.pars)[j],sep=""),type="l",lwd=1,lty=1)
}


#ACF plot
burnin<-1000
par(mfrow=c(p,1))
for(j in 1:p){
  acf(mcmc3.pars[(burnin+1):S,j],main=paste("ACF for ",colnames(mcmc3.pars)[j],sep=""),lag.max=100)
}


thin<-25
par(mfrow=c(p,1))
for(j in 1:p){# after thinning
  acf(mcmc3.pars[seq(burnin+1,S,by=thin),j],main=paste("ACF for ",colnames(mcmc3.pars)[j],sep=""),lag.max=100)
}

length(mcmc3.pars[seq(burnin+1,S,by=thin),1])
for(j in 1:p){
  print(effectiveSize(mcmc3.pars[seq(burnin+1,S,by=thin),j]))
}


mcmc3.pars.pMean<-apply(mcmc3.pars[seq(burnin+1,S,by=thin),],2,mean) # posterior mean of parameters
mcmc3.pars.pMean
mcmc3.pars.pSD<-apply(mcmc3.pars[seq(burnin+1,S,by=thin),],2,sd) # posterior mean of parameters
mcmc3.pars.pSD




