library(MASS)
library(msm)
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
data<-na.omit(cbind(y,X)) # remove missing data
y<-data[,1]
X<-data[,-1]

ranks<-match(y,sort(unique(y))) 
uranks<-sort(unique(ranks))
K<-length(uranks)

n<-dim(X)[1] ; p<-dim(X)[2]
XX<-t(X)%*%X  

#### Ordinal probit regression
## setup
beta.0<-rep(0,p) # prior mean of beta
Sigma.0<-n*solve(XX) # prior covariance matrix of beta
iSigma.0<-solve(Sigma.0)
mu.0<-rep(0,K-1) # prior mean of g
gamma.0<-rep(10,K-1) #prior SD of g

beta<-rep(0,p) #initial values
z<-qnorm(rank(y,ties.method="random")/(n+1)) #initial values


g<-rep(NA,length(uranks)-1)


## MCMC
S<-26000
beta.MCMC<-matrix(NA,S,p) 
z.MCMC<-matrix(NA,S,n) 
g.MCMC<-matrix(NA,S,length(uranks)-1)

ac<-0

set.seed(0)
for(s in 1:S) 
{

  #update g 
  for(k in 1:(K-1)) 
  {
    a<-max(z[y==k])
    b<-min(z[y==(k+1)])
    g[k]<- rtnorm(1,mu.0[k],gamma.0[k],a,b)
  }

  #update beta
  Sigma.n<-solve(XX+iSigma.0) # posterior covariance matrix of beta
  beta.n<-Sigma.n%*%(t(X)%*%z+iSigma.0%*%beta.0) # posterior mean of beta
  beta<-mvrnorm(1,beta.n,Sigma.n)
  
  
  #update z
  ez<-X%*%beta 
  a<-c(-Inf,g)[ match( y-1, 0:K) ]
  b<-c(g,Inf)[y]  
  z<-rtnorm(n, ez, 1, a, b) 
  
  
  #help mixing, based on the code from http://www2.stat.duke.edu/~pdh10/FCBS/Replication/chapter12.R
  c<-rnorm(1,0,n^(-1/3))  
  zp<-z+c ; gp<-g+c
  lhr<-  sum(dnorm(zp,ez,1,log=T) - dnorm(z,ez,1,log=T) ) + 
    sum(dnorm(gp,mu.0,gamma.0,log=T) - dnorm(g,mu.0,gamma.0,log=T) )
  if(log(runif(1))<lhr) { z<-zp ; g<-gp ; ac<-ac+1 }

  
  beta.MCMC[s,]<-beta
  z.MCMC[s,]<-z
  g.MCMC[s,]<-g
  
  if(s%%1000==0){
    print(s)
  }
} 


# Traceplot
par(mfrow=c(p,1))
for(j in 1:p){# plot the first 2000 iterations
  plot(beta.MCMC[1:2000,j],ylab=paste("beta_",j,sep=""),xlab="Iteration",main=paste("Trace plot for beta_",j,sep=""),type="l",lwd=1,lty=1)
}


#ACF plot
burnin<-1000
par(mfrow=c(p,1))
for(j in 1:p){
  acf(beta.MCMC[(burnin+1):S,j],main=paste("ACF for beta_",j,sep=""),lag.max=100)
}


thin<-25
par(mfrow=c(p,1))
for(j in 1:p){# after thinning
  acf(beta.MCMC[seq(burnin+1,S,by=thin),j],main=paste("ACF for beta_",j,sep=""),lag.max=100)
}

length(beta.MCMC[seq(burnin+1,S,by=thin),1])
for(j in 1:p){
print(effectiveSize(beta.MCMC[seq(burnin+1,S,by=thin),j]))
}


beta.pMean<-apply(beta.MCMC[seq(burnin+1,S,by=thin),],2,mean) # posterior mean of beta
beta.pMean
beta.pSD<-apply(beta.MCMC[seq(burnin+1,S,by=thin),],2,sd) # posterior mean of beta
beta.pSD

g.pMean<-apply(g.MCMC[seq(burnin+1,S,by=thin),],2,mean)
g.pMean
g.pSD<-apply(g.MCMC[seq(burnin+1,S,by=thin),],2,sd)
g.pSD

#### Figure 12.2 in the textbok
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,2))
plot(X[,1]+.25*(X[,2]),z.MCMC[S,],
 pch=15+X[,2],col=c("gray","black")[X[,2]+1],
 xlab="number of children",ylab="z", ylim=range(c(-2.5,4,z.MCMC[S,])),
    xlim=c(0,9))

abline(0,beta.pMean[1],lwd=2 ,col="gray")
abline(beta.pMean[2],beta.pMean[1]+beta.pMean[3],col="black",lwd=2 )
legend(5,4,legend=c("PDEG=0","PDEG=1"),pch=c(15,16),col=c("gray","black"))


plot(density(beta.MCMC[seq(burnin+1,S,by=thin),3]),lwd=2,xlim=c(-.5,.5),main="",
    xlab=expression(beta[3]),ylab="density")
sd<-sqrt(  solve(t(X)%*%X/n)[3,3] )
x<-seq(-.7,.7,length=100)
lines(x,dnorm(x,0,sd),lwd=2,col="gray")
legend(-.5,6,legend=c("prior","posterior"),lwd=c(2,2),col=c("gray","black"),bty="n")












