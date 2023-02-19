#####  Fitting a hierarchical linear regression model
#####  1st stage: y_j| beta_j, X_j, sigma2 ~ N_nj(beta^prime_j X_j, sigma2 I_nj) for j=1,..,m  
#####  2nd stage: beta_j | theta, Sigma ~ N_p(theta,Sigma)  for j=1,..,m
#####  3rd stage (priors):    theta ~ N_p(mu_0, Lambda_0) 
#####                         Sigma ~ InverseWishart(eta_0,S_0)
#####                         sigma2 ~ InverseGamma(nu0/2,nu0*sigma2_0/2)


library(MASS)
library(coda)
library(MCMCpack)

## Data

math.data <- read.table("math_scores_ses.txt",header=TRUE)
math.data[1:3,]

school <- math.data$school_id
m <- length(unique(school))
y <- math.data$math_score
n <- length(y)
ses <- math.data$SES
p <- 2
X <- matrix(cbind(rep(1,n),ses),nrow=n,ncol=p)

## determining the number of observations for each group
n.j <- as.numeric(table(school))


## OLS estimates of regression coefficients
beta.ols <- matrix(0,nrow=m,ncol=p)
var.school <- rep(0,m)
for(j in 1:m){
	y.j <- y[which(school==unique(school)[j])]
	ses.j <- ses[which(school==unique(school)[j])]
	reg.j <- lm(y.j~ses.j)
	beta.ols[j,] <- as.numeric(reg.j$coeff)
	var.school[j] <- var(y.j)
}

plot(seq(-3,3,length=100),beta.ols[1,1]+beta.ols[1,2]*seq(-3,3,length=100),xlab="SES",ylab="Math score",
col="grey",type="l",xlim=c(-3,3),ylim=c(20,85),lwd=2,cex=2,cex.axis=2,cex.lab=1.5)
for(j in 2:m){
 lines(seq(-3,3,length=100),beta.ols[j,1]+beta.ols[j,2]*seq(-3,3,length=100),col="grey",lwd=2)
}
lines(seq(-3,3,length=100),mean(beta.ols[,1])+mean(beta.ols[,2])*seq(-3,3,length=100),col="black",lwd=3)


plot(n.j,beta.ols[,1],xlab="Sample size n_j",ylab="OLS estimate of \n Intercept",type="p",pch=20,col="black",cex=2,cex.axis=2,cex.lab=1.5)
abline(h=mean(beta.ols[,1]),col="red",lwd=2)

plot(n.j,beta.ols[,2],xlab="Sample size",ylab="OLS estimate of \n Slope",type="p",pch=20,col="black",cex=2,cex.axis=2,cex.lab=1.5)
abline(h=mean(beta.ols[,2]),col="red",lwd=2)


## Parameters of prior distributions
## prior on theta
mu0 <- colMeans(beta.ols)
Lambda0 <- cov(beta.ols)

## prior on Sigma
eta0 <- p+2
S0 <- Lambda0
## prior on sigma^2
nu0 <- 1
sigma2.0 <- mean(var.school)





## Functions to sample from full conditionals

sample.beta<- function(y,X,n.j,p,m,group,theta,Sigma,sigma2){

	new.beta <- matrix(0,m,p)
	for(j in 1:m){
		
		# subsetting data
		y.j <- y[group==unique(group)[j]]
		X.j <- X[group==unique(group)[j],]
		
		post.var <- solve(solve(Sigma)+((t(X.j)%*%X.j)/sigma2))
		post.mean <- post.var%*%((solve(Sigma)%*%matrix(theta,nrow=p,ncol=1)) + ((t(X.j)%*%matrix(y.j,nrow=n.j[j],ncol=1))/sigma2))
		new.beta[j,] <- mvrnorm(1,post.mean,post.var)
	}	
	
	return(new.beta)
}



sample.theta <- function(beta,m,p,Sigma,mu0,Lambda0){

	beta.bar <- as.numeric(apply(beta,2,mean))
	
	post.var <- solve(solve(Lambda0)+m*solve(Sigma))
	post.mean <- post.var%*%((solve(Lambda0)%*%matrix(mu0,nrow=p,ncol=1))+m*(solve(Sigma)%*%matrix(beta.bar,nrow=p,ncol=1)))
	new.theta <- mvrnorm(1,post.mean,post.var)
	return(new.theta)
}


sample.Sigma <- function(beta,theta,m,p,eta0,S0){
	
	eta.n <- eta0+m
	S.theta <- matrix(0,nrow=p,ncol=p)
	for(j in 1:m){
		S.theta <- S.theta+(beta[j,]-theta)%*%t(beta[j,]-theta)
	}
	S.theta <- S0+S.theta
	new.Sigma <- solve(rwish(eta.n,solve(S.theta)))
	return(new.Sigma)
}


sample.sigma2 <- function(y,X,beta,group,m,n,n.j,nu0,sigma2.0){
	
	nu.n <- nu0+n
	sse <- 0
	for(j in 1:m){
		y.j <- y[group==unique(group)[j]]
		X.j <- X[group==unique(group)[j],]
		sse <- sse + (t(matrix(y.j,nrow=n.j[j],ncol=1) -X.j%*%matrix(beta[j,],nrow=p,ncol=1))%*%(matrix(y.j,nrow=n.j[j],ncol=1) -X.j%*%matrix(beta[j,],nrow=p,ncol=1)))
	}
	sigma2.n <- (1/nu.n)*((nu0*sigma2.0)+sse)
	new.sigma2 <- 1/rgamma(1,nu.n/2,(nu.n*sigma2.n)/2)
	return(new.sigma2)
}




### Running the Metropolis algorithm

S <- 10000


## Initial values 
init.beta <- beta.ols
init.sigma2 <- mean(var.school)
init.Sigma <- Lambda0
init.theta <- apply(beta.ols,2,mean)


## We store the MCMC samples for beta in an array with m rows, 2 columns and S layers, one for each MCMC iteration 
beta.MCMC <- array(0,dim=c(m,p,S))
## We store the MCMC samples for sigma2 in a S times 1 vector  
sigma2.MCMC <- rep(0,S)
## We store the MCMC samples for theta in a S times p matrix
theta.MCMC <- matrix(0,S,p)
## We store the MCMC samples for Sigma in a matrix with S rows and p+(p*(p-1)/2) (here p=2)
Sigma.MCMC <- matrix(0,S,3)

## We store the predicted beta1 and beta2 for a new school yet to be sampled in a S times p matrix
pred.beta <- matrix(0,S,p)

set.seed(0)
for(k in 1:S){
	
	if(k==1){
		beta <- init.beta
		sigma2 <- init.sigma2
		Sigma <- init.Sigma
		theta <- init.theta
	}
	

	new.theta <- sample.theta(beta,m,p,Sigma,mu0,Lambda0)
	new.Sigma <- sample.Sigma(beta,new.theta,m,p,eta0,S0)
	new.sigma2 <- sample.sigma2(y,X,beta,school,m,n,n.j,nu0,sigma2.0)
	new.beta <- sample.beta(y,X,n.j,p,m,school,new.theta,new.Sigma,new.sigma2)
	

	theta <- new.theta
	Sigma <- new.Sigma
	sigma2 <- new.sigma2
	beta <- new.beta
	
	theta.MCMC[k,] <- theta
	Sigma.MCMC[k,] <- c(Sigma[1,1],Sigma[2,2],Sigma[1,2])
	sigma2.MCMC[k] <- sigma2
	beta.MCMC[,,k] <- beta
	pred.beta[k,] <- mvrnorm(1,theta,Sigma)
	
	print(k)
}



## Effective sample size
effectiveSize(sigma2.MCMC)
effectiveSize(theta.MCMC[,1])
effectiveSize(theta.MCMC[,2])
effectiveSize(Sigma.MCMC[,1])
effectiveSize(Sigma.MCMC[,2])
effectiveSize(Sigma.MCMC[,3])


## Trace plots
par(mfrow=c(2,2))
for(j in 1:2){
	plot(beta.MCMC[j,1,],ylab=paste("Beta_1,",j,sep=""),xlab="Iteration",main=paste("Trace plot for beta_1 \n School=",j,sep=""),type="l",lwd=1,lty=1)
	plot(beta.MCMC[j,2,],ylab=paste("Beta_2,",j,sep=""),xlab="Iteration",main=paste("Trace plot for beta_2 \n School=",j,sep=""),type="l",lwd=1,lty=1)
}

par(mfrow=c(1,3))
plot(Sigma.MCMC[,1],ylab="Sigma_1,1",xlab="Iteration",main="Trace plot for Sigma_1,1",type="l",lwd=1,lty=1)
plot(Sigma.MCMC[,2],ylab="Sigma_2,2",xlab="Iteration",main="Trace plot for Sigma_2,2",type="l",lwd=1,lty=1)
plot(Sigma.MCMC[,3],ylab="Sigma_1,2",xlab="Iteration",main="Trace plot for Sigma_1,2",type="l",lwd=1,lty=1)

par(mfrow=c(1,3))
plot(sigma2.MCMC,ylab="sigma^2",xlab="Iteration",main="Trace plot for sigma^2",type="l",lwd=1,lty=1)
plot(theta.MCMC[,1],ylab="Theta_1",xlab="Iteration",main="Trace plot for Theta_1",type="l",lwd=1,lty=1)
plot(theta.MCMC[,2],ylab="Theta_2",xlab="Iteration",main="Trace plot for Theta_2",type="l",lwd=1,lty=1)


## Marginal posterior densities

thinning <- 10

par(mfrow=c(1,1))
plot(density(theta.MCMC[seq(1,S,by=thinning),1]),xlab="Theta_1",ylab="Density",
main="Marginal posterior density \n for Theta_1",col="black",lwd=2,xlim=c(40,60))
lines(density(pred.beta[seq(1,S,by=thinning),1]),col="blue",lwd=2)


par(mfrow=c(1,1))
plot(density(theta.MCMC[seq(1,S,by=thinning),2]),xlab="Theta_2",ylab="Density",
main="Marginal posterior density \n for Theta_2",col="black",lwd=2,xlim=c(-8,8))
lines(density(pred.beta[seq(1,S,by=thinning),2]),col="blue",lwd=2)



par(mfrow=c(1,1))
plot(density(sigma2.MCMC[seq(1,S,by=thinning)]),xlab="sigma^2",ylab="Density",main="Marginal posterior density \n for sigma^2",col="black",lwd=2)


par(mfrow=c(1,3))
plot(density(Sigma.MCMC[seq(1,S,by=thinning),1]),xlab="Sigma_1,1",ylab="Density",main="Marginal posterior density \n for Sigma_1,1",col="black",lwd=2)
plot(density(Sigma.MCMC[seq(1,S,by=thinning),2]),xlab="Sigma_2,2",ylab="Density",main="Marginal posterior density \n for Sigma_2,2",col="black",lwd=2)
plot(density(Sigma.MCMC[seq(1,S,by=thinning),3]),xlab="Sigma_1,2",ylab="Density",main="Marginal posterior density \n for Sigma_1,2",col="black",lwd=2)



## Posterior summaries: burn-in=1000
post.mean.beta <- matrix(0,m,2)
post.sd.beta <- matrix(0,m,2)
post.l95.beta <- matrix(0,m,2)
post.u95.beta <- matrix(0,m,2)

for(j in 1:m){
	post.mean.beta[j,1] <- mean(beta.MCMC[j,1,seq(1000,S,by=thinning)])
	post.mean.beta[j,2] <- mean(beta.MCMC[j,2,seq(1000,S,by=thinning)])
	post.sd.beta[j,1] <- sd(beta.MCMC[j,1,seq(1000,S,by=thinning)])
	post.sd.beta[j,2] <- sd(beta.MCMC[j,2,seq(1000,S,by=thinning)])
	post.l95.beta[j,1] <- sd(beta.MCMC[j,1,seq(1000,S,by=thinning)],0.025)
	post.l95.beta[j,2] <- sd(beta.MCMC[j,2,seq(1000,S,by=thinning)],0.025)
	post.u95.beta[j,1] <- quantile(beta.MCMC[j,1,seq(1000,S,by=thinning)],0.975)
	post.u95.beta[j,2] <- quantile(beta.MCMC[j,2,seq(1000,S,by=thinning)],0.975)
}

post.mean.sigma2 <- mean(sigma2.MCMC[seq(1000,S,by=thinning)])
post.sd.sigma2 <- sd(sigma2.MCMC[seq(1000,S,by=thinning)])
post.95ci.sigma2 <- quantile(sigma2.MCMC[seq(1000,S,by=thinning)],c(0.025,0.975))

post.mean.theta <- apply(theta.MCMC[seq(1000,S,by=thinning),],2,mean)
post.sd.theta <- apply(theta.MCMC[seq(1000,S,by=thinning),],2,sd)
post.95ci.theta <- apply(theta.MCMC[seq(1000,S,by=thinning),],2,quantile,c(0.025,0.975))


post.mean.Sigma <- apply(Sigma.MCMC[seq(1000,S,by=thinning),],2,mean)
post.sd.Sigma <- apply(Sigma.MCMC[seq(1000,S,by=thinning),],2,sd)
post.95ci.Sigma <- apply(Sigma.MCMC[seq(1000,S,by=thinning),],2,quantile,c(0.025,0.975))




