
library(MASS)
library(coda)

#### Gibbs sampling for the linear regression model with semi-conjugate prior p(beta, sigma^2)=N(beta; beta0,Sigma0)*InverseGamma(nu0/2,nu0*(sigma02)/2)

# Data on the maximal oxygen uptake example
y <- c(-0.87, -10.74, -3.27, -1.97, 7.5, -7.25, 17.05, 4.96, 
10.4, 11.05, 0.26, 2.51)
program <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1)
age <- c(23, 22, 22, 25, 27, 20, 31, 
23, 27, 28, 22, 24)
program.age <- c(0, 0, 0, 0, 0, 0, 31, 23, 27, 28, 22, 24)

## Setting up the matrices
n <- length(y)
p <- 4
X <- matrix(cbind(rep(1,n),program,age,program.age),nrow=n,ncol=p)

Xprime.X <- t(X)%*%X
Xprime.y <- t(X)%*%matrix(y,nrow=n,ncol=1)

## OLS estimates
beta.ols <- solve(Xprime.X)%*%Xprime.y
sigma2.ols <- (1/(n-p))*t(y-X%*%beta.ols)%*%(y-X%*%beta.ols)



## Prior parameters
## prior parameters for beta
beta0 <- c(0,0,0,0)
Sigma0 <- matrix(0,nrow=p,ncol=p)
Sigma0[1,1] <- 150^2
Sigma0[2,2] <- 150^2
Sigma0[3,3] <- 6^2
Sigma0[4,4] <- 6^2

## prior precision matrix
iSigma0 <- solve(Sigma0)
## prior precision times prior mean
iSigma0beta0 <- iSigma0%*%matrix(beta0,nrow=p,ncol=1)

## prior parameters for sigma2
nu0 <- 1
sigma2.0 <- 40

nu.n <- nu0+n

## Covariates for a new observation
m <- 2
Xtilde <- matrix(c(1,0,30,0,1,1,30,30),nrow=m,ncol=p,byrow=T)

### prediction with OLS
pred.ols <- Xtilde%*%beta.ols

### functions to sample from full conditionals
sample.beta <- function(y,sigma2,Xprime.X,Xprime.y,iSigma0,iSigma0beta0){
	
	post.var <- solve(iSigma0+(Xprime.X/sigma2))
	post.mean <- post.var%*%(iSigma0beta0+(Xprime.y/sigma2))
	new.beta <- mvrnorm(1,post.mean,post.var)
	return(new.beta)
}
sample.sigma2 <- function(y,X,beta,nu.n,nu0,sigma2.0,p,n){
	
	ss.residuals <- t(y-X%*%matrix(beta,nrow=p,ncol=1))%*%(y-X%*%matrix(beta,nrow=p,ncol=1))
	sigma2.n <- 1/(nu.n)*((nu0*sigma2.0)+ss.residuals)
	new.sigma2 <- 1/rgamma(1,nu.n/2,(nu.n*sigma2.n)/2)
	return(new.sigma2)
}

## Running the Gibbs sampling algorithm
S <- 10000
output.MCMC <- matrix(0,S,p+1)
pred.MCMC <- matrix(0,S,m)

## Initial values
init.beta <- c(0,0,0,0)
init.sigma2 <- var(y)

set.seed(0)
for(i in 1:S){
	if(i==1){
		output.MCMC[i,1] <- init.sigma2
		output.MCMC[i,2:(p+1)] <- init.beta
		pred.MCMC[i,1:m] <- mvrnorm(1,Xtilde%*%output.MCMC[i,2:(p+1)],output.MCMC[i,1]*diag(1,m))
	}
	if(i>1){
		output.MCMC[i,1] <- sample.sigma2(y,X,output.MCMC[(i-1),(2:(p+1))],nu.n,nu0,sigma2.0,p,n)
		output.MCMC[i,2:(p+1)] <- sample.beta(y,output.MCMC[i,1],Xprime.X,Xprime.y,iSigma0,iSigma0beta0)
		pred.MCMC[i,1:m] <- mvrnorm(1,Xtilde%*%output.MCMC[i,2:(p+1)],output.MCMC[i,1]*diag(1,m))
	}
	print(i)
}

### Examining trace plots, autocorrelation functions and effective sample size
par(mfrow=c(2,2))
plot(output.MCMC[,2],xlab="Iteration number",ylab="Beta_1",type="l",col="black")
plot(output.MCMC[,3],xlab="Iteration number",ylab="Beta_2",type="l",col="black")
plot(output.MCMC[,4],xlab="Iteration number",ylab="Beta_3",type="l",col="black")
plot(output.MCMC[,5],xlab="Iteration number",ylab="Beta_4",type="l",col="black")

par(mfrow=c(1,1))
plot(output.MCMC[,1],xlab="Iteration number",ylab="Sigma2",type="l",col="black")

par(mfrow=c(2,2))
acf(output.MCMC[,2],main="Beta_1")
acf(output.MCMC[,3],main="Beta_2")
acf(output.MCMC[,4],main="Beta_3")
acf(output.MCMC[,5],main="Beta_4")

par(mfrow=c(1,1))
acf(output.MCMC[,1],main="Sigma2")

effectiveSize(output.MCMC[,1])
effectiveSize(output.MCMC[,2])
effectiveSize(output.MCMC[,3])
effectiveSize(output.MCMC[,4])
effectiveSize(output.MCMC[,5])


### Marginal posterior densities for the regression parameters
burnin <- 1000

#pdf("beta_oxygen_semiconjugate_prior_pdfs.pdf",width=6,height=6)
par(mfrow=c(2,2))
plot(density(output.MCMC[(burnin+1):S,2]),col="black",main="Beta_1")
abline(v=mean(output.MCMC[(burnin+1):S,2]),col="red",lwd=2,lty=1)
abline(v=beta.ols[1],col="blue",lwd=2,lty=1)
abline(v=beta0[1],col="grey",lwd=2,lty=1)
abline(v=quantile(output.MCMC[(burnin+1):S,2],0.025),lty=3,lwd=2,col="red")
abline(v=quantile(output.MCMC[(burnin+1):S,2],0.975),lty=3,lwd=2,col="red")

plot(density(output.MCMC[(burnin+1):S,3]),col="black",main="Beta_2")
abline(v=mean(output.MCMC[(burnin+1):S,3]),col="red",lwd=2,lty=1)
abline(v=beta.ols[2],col="blue",lwd=2,lty=1)
abline(v=beta0[2],col="grey",lwd=2,lty=1)
abline(v=quantile(output.MCMC[(burnin+1):S,3],0.025),lty=3,lwd=2,col="red")
abline(v=quantile(output.MCMC[(burnin+1):S,3],0.975),lty=3,lwd=2,col="red")

plot(density(output.MCMC[(burnin+1):S,4]),col="black",main="Beta_3")
abline(v=mean(output.MCMC[(burnin+1):S,4]),col="red",lwd=2,lty=1)
abline(v=beta.ols[3],col="blue",lwd=2,lty=1)
abline(v=beta0[3],col="grey",lwd=2,lty=1)
abline(v=quantile(output.MCMC[(burnin+1):S,4],0.025),lty=3,lwd=2,col="red")
abline(v=quantile(output.MCMC[(burnin+1):S,4],0.975),lty=3,lwd=2,col="red")

plot(density(output.MCMC[(burnin+1):S,5]),col="black",main="Beta_4")
abline(v=mean(output.MCMC[(burnin+1):S,5]),col="red",lwd=2,lty=1)
abline(v=beta.ols[4],col="blue",lwd=2,lty=1)
abline(v=beta0[4],col="grey",lwd=2,lty=1)
abline(v=quantile(output.MCMC[(burnin+1):S,5],0.025),lty=3,lwd=2,col="red")
abline(v=quantile(output.MCMC[(burnin+1):S,5],0.975),lty=3,lwd=2,col="red")
#dev.off()


#pdf("sigma2_oxygen_semiconjugate_prior_pdfs.pdf",width=6,height=6)
par(mfrow=c(1,1))
plot(density(output.MCMC[(burnin+1):S,1]),col="black",main="Sigma^2")
abline(v=mean(output.MCMC[(burnin+1):S,1]),col="red",lwd=2,lty=1)
abline(v=sigma2.ols,col="blue",lwd=2,lty=1)
abline(v=sigma2.0,col="grey",lwd=2,lty=1)
abline(v=quantile(output.MCMC[(burnin+1):S,1],0.025),lty=3,lwd=2,col="red")
abline(v=quantile(output.MCMC[(burnin+1):S,1],0.975),lty=3,lwd=2,col="red")
#dev.off()

#pdf("pred_oxygen_semiconjugate_prior_pdfs.pdf",width=6,height=6)
par(mfrow=c(2,1))
plot(density(pred.MCMC[(burnin+1):S,1]),col="black",main="Prediction for new individual 1")
abline(v=mean(pred.MCMC[(burnin+1):S,1]),col="red",lwd=2,lty=1)
abline(v=pred.ols[1],col="blue",lwd=2,lty=1)
abline(v=quantile(pred.MCMC[(burnin+1):S,1],0.025),lty=3,lwd=2,col="red")
abline(v=quantile(pred.MCMC[(burnin+1):S,1],0.975),lty=3,lwd=2,col="red")

plot(density(pred.MCMC[(burnin+1):S,2]),col="black",main="Prediction for new individual 2")
abline(v=mean(pred.MCMC[(burnin+1):S,2]),col="red",lwd=2,lty=1)
abline(v=pred.ols[2],col="blue",lwd=2,lty=1)
abline(v=quantile(pred.MCMC[(burnin+1):S,2],0.025),lty=3,lwd=2,col="red")
abline(v=quantile(pred.MCMC[(burnin+1):S,2],0.975),lty=3,lwd=2,col="red")
#dev.off()

### Posterior summaries for the parameters

# Posterior mean
apply(output.MCMC[(burnin+1):S,],2,mean)

# Posterior median
apply(output.MCMC[(burnin+1):S,],2,median)

# Posterior standard deviation
apply(output.MCMC[(burnin+1):S,],2,sd)


### Posterior summaries for the prediction

# Posterior mean or posterior predictive mean
apply(pred.MCMC[(burnin+1):S,1:2],2,mean)

# Posterior median
apply(pred.MCMC[(burnin+1):S,1:2],2,median)

# Posterior standard deviation
apply(pred.MCMC[(burnin+1):S,1:2],2,sd)



###  Model checking

res.MCMC <- matrix(0,S,n)
PPO.MCMC <- matrix(0,S,n)

for(i in 1:S){
	res.MCMC[i,1:n] <- y-X%*%output.MCMC[i,2:(p+1)]
	PPO.MCMC[i,1:n] <- 1-pnorm(y,X%*%output.MCMC[i,2:(p+1)],sqrt(output.MCMC[i,1]))
}

### Probability that the i-th observation is a residual

par(mfrow=c(4,3))
for(i in 1:12){
	plot(density(res.MCMC[(burnin+1):S,i]),main=paste("residual_",i,sep=""),col="black")	
}

## Posterior mean for the standard deviation sigma
sigma.mean <- mean(sqrt(output.MCMC[(burnin+1):S,1]))

## Posterior probability that a residual is larger than 3 times sigma
apply(abs(res.MCMC[(burnin+1):S,]) > 3*sigma.mean,2,mean)

apply(PPO.MCMC[(burnin+1):S,],2,mean)
# Note that here PPO_i is computed as the posterior prob > y_i instead of the posterior density at y_i, 
# so if PPO_i is beyond [0.025,0.975] (or other user-prespecified interval), 
# then y_i can be viewed as an outlier.
