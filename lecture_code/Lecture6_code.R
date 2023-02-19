### Code to perform Gibbs sampling for the multivariate normal data
### sampling model: y1,..,yn ~ N_2(theta,Sigma)
### priors:         p(theta,Sigma)=p(theta)*p(Sigma)
###                 p(theta)=N_2(mu0, Lambda0)
###                 p(Sigma)=Inverse Wishart(nu0,Phi0)

library(MCMCpack)
library(MASS)

## Data

y <- matrix(c(59, 43, 34, 32, 42, 38, 55, 67, 64, 45, 49, 72, 34, 
70, 34, 50, 41, 52, 60, 34, 28, 35, 77, 39, 46, 26, 38, 43, 68, 
86, 77, 60, 50, 59, 38, 48, 55, 58, 54, 60, 75, 47, 48, 33),nrow=22,ncol=2,byrow=FALSE)


## number of observations
n <- dim(y)[1]

## ybar vector
ybar <- apply(y,2,mean)


## matrix with sample variances and covariances
S2 <- cov(y)


## parameters of prior distribution
mu0 <- c(50,50)
Lambda0 <- matrix(c(625,312.5,312.5,625),nrow=2,ncol=2,byrow=T)

nu0 <- 4
Phi0 <- Lambda0


### functions to do the sampling from the full conditionals

## To sample from the multivariate normal distribution we will need to use the package MASS 
## that contains the mvrnorm function to sample from a multivariate normal distribution
## and the R package MCMCpack that contains the rwish function 
## to sample from a Wishart distribution


sample.theta <- function(ybar,mu0,Lambda0,n,Sigma){

	# prior precision matrix
	iLambda0 <- solve(Lambda0)
	
	# data precision matrix
	iSigma <- solve(Sigma)
	
	# posterior precision matrix
	post.prec <- n*iSigma + iLambda0
	
	# posterior variance matrix
	post.var <- solve(post.prec)
	
	# posterior mean
	post.mean <- post.var%*%(iLambda0%*%mu0 + n*iSigma%*%ybar)
	
	new.theta <- mvrnorm(1,post.mean,post.var)
	return(new.theta)
}


sample.Sigma <- function(y,theta,n,nu0,Phi0){
	
	# posterior degrees of freedom
	nu.n <- n+nu0
	
	
	# residual sum of squares
	rss <- (t(y)-theta)%*%t(t(y)-theta)
	# Phi.n
	Phi.n <- Phi0+rss
	
	new.iSigma <- rwish(nu.n,solve(Phi.n))
	new.Sigma <- solve(new.iSigma)
	return(new.Sigma)
}



### Gibbs sampling algorithm

S <- 10000
theta.MC <- matrix(0,S,2)
Sigma.MC <- array(0,dim=c(2,2,S))

set.seed(0)
# set the random seed for reproducible work! 
# Note that I used a different seed for the lecture note. Compare the results.
for(i in 1:B){
	if(i==1){
		theta.MC[i,] <- ybar
		Sigma.MC[,,i] <- S2
	}
	if(i>1){
		theta.MC[i,] <- sample.theta(ybar,mu0,Lambda0,n,Sigma.MC[,,(i-1)])	
		Sigma.MC[,,i] <- sample.Sigma(y,theta.MC[i,],n,nu0,Phi0)
	}
	
    print(i)
		
}


###  Trace plot

par(mfrow=c(1,2))
plot(theta.MC[,1],xlab="Iteration number",ylab="Theta_1",type="l")
plot(theta.MC[,2],xlab="Iteration number",ylab="Theta_2",type="l")


par(mfrow=c(1,3))
plot(Sigma.MC[1,1,],xlab="Iteration number",ylab="sigma^2_1",type="l")
plot(Sigma.MC[1,2,],xlab="Iteration number",ylab="sigma_12",type="l")
plot(Sigma.MC[2,2,],xlab="Iteration number",ylab="sigma^2_2",type="l")


###  ACF plot

par(mfrow=c(1,2))
acf(theta.MC[,1],main="Theta_1")
acf(theta.MC[,2],main="Theta_2")


par(mfrow=c(1,3))
acf(Sigma.MC[1,1,],main="Sigma^2_1")
acf(Sigma.MC[1,2,],main="Sigma^2_2")
acf(Sigma.MC[2,2,],main="Sigma_12")


## Posterior summaries

## difference between the after and before test scores

burnin <- 5000
quantile(theta.MC[((burnin+1):S),2]-theta.MC[((burnin+1):S),1],c(0.025,0.975))
mean(theta.MC[((burnin+1):S),2]-theta.MC[((burnin+1):S),1])
median(theta.MC[((burnin+1):S),2]-theta.MC[((burnin+1):S),1])
sd(theta.MC[((burnin+1):S),2]-theta.MC[((burnin+1):S),1])

## correlation between the after and before test scores

quantile(Sigma.MC[1,2,((burnin+1):S)]/sqrt(Sigma.MC[1,1,((burnin+1):S)]*Sigma.MC[2,2,((burnin+1):S)]),c(0.025,0.975))
mean(Sigma.MC[1,2,((burnin+1):S)]/sqrt(Sigma.MC[1,1,((burnin+1):S)]*Sigma.MC[2,2,((burnin+1):S)]))
median(Sigma.MC[1,2,((burnin+1):S)]/sqrt(Sigma.MC[1,1,((burnin+1):S)]*Sigma.MC[2,2,((burnin+1):S)]))
sd(Sigma.MC[1,2,((burnin+1):S)]/sqrt(Sigma.MC[1,1,((burnin+1):S)]*Sigma.MC[2,2,((burnin+1):S)]))




									 
									 
									 
									 












