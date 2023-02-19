####  Hierarchical generalized linear model:
####  Data: y_ij, j=1,..,m, i=1,..,n_j
####  1st stage model: y_ij | theta_ij ~ Poisson(theta_ij) for j=1,..,m; i=1,..,n_j
####                   theta_ij = exp(beta_j * x_i)    
####                   x_i=(1 i/20 (i/20)^2  (i/20)^3 (i/20)^4)
####  2nd stage model: beta_j | mu, Sigma ~ N_5(mu, Sigma) for j=1,..,m
####  Priors:          mu ~ N_5(mu_0,Lambda_0)
####                   Lambda_0 ~ InverseWishart(eta_0, S^{-1}_0)

####  Gibbs sampling algorithm with Metropolis step
####  (for forms of full conditionals see Lecture 21)

library(coda)
library(MCMCpack)
library(MASS)
library(mvtnorm)

#############
####  Data
#############

tumor.data <- read.table("Tumor_counts.txt",header=FALSE)
Y <- as.matrix(tumor.data)

## Each row is a mouse
dim(tumor.data)

## Number of observations per mouse
n.j <- dim(tumor.data)[2]
## Number of mice
m <- dim(tumor.data)[1]


## Plotting the data
par(mfrow=c(1,1))

plot(1:20,Y[1,],type="l",xlab="Location",ylab="Number of tumors",lwd=2,lty=1,col="grey",ylim=c(0,18))
for(j in 2:m){
	lines(1:20,Y[j,],type="l",lwd=2,lty=1,col="grey")
}
lines(1:20,apply(Y,2,mean),col="black",lwd=2,lty=1,type="l")


## Creating the X matrix: each tumor number is a function of the tumor location
## and we model the relationship using a fourth degree polynomial
p <- 5
X <- matrix(0,nrow=n.j,ncol=5)
for(i in 1:n.j){
	X[i,] <- c(1,(i/20),(i/20)^2,(i/20)^3,(i/20)^4)
}


########################
### Prior distributions
########################

## Fitting a separate linear regression on the log average number of tumor per mouse
beta.ols <- matrix(0,nrow=m,ncol=5)
for(j in 1:m){
	log.data.mouse <- rep(0,20)
	for(i in 1:20){
		log.data.mouse[i] <- log(Y[j,i]+(i/20))	
	}
	beta.ols[j,] <- solve(t(X)%*%X)%*%(t(X)%*%matrix(log.data.mouse,nrow=20,ncol=1))	
}

mu0 <- apply(beta.ols,2,mean) # or colMeans(beta.ols)
Lambda0 <- cov(beta.ols)

eta0 <- 7
S0 <- Lambda0


########################
### Sampling functions
########################

sample.beta.j <- function(current.beta.j,delta,X,Yj,mu,Sigma,nj){
	
	# proposing a candidate value
	beta.j.star <- mvrnorm(1,current.beta.j,delta*Sigma)
	accept <- 0
	
	# evaluating ratio on the log scale
	log.num <- 0
	log.den <- 0
	for(i in 1:nj){
		theta.ij.star <- exp(sum(beta.j.star*X[i,]))
		theta.ij <- exp(sum(current.beta.j*X[i,]))
		log.num <- log.num+dpois(Yj[i],theta.ij.star,log=TRUE)
		log.den <- log.den+dpois(Yj[i],theta.ij,log=TRUE)
	}
	log.num <- dmvnorm(beta.j.star,mu,Sigma,log=TRUE)+log.num
	log.den <- dmvnorm(current.beta.j,mu,Sigma,log=TRUE)+log.den
	
	
	r <- exp(log.num-log.den)
	u <- runif(1,0,1)
	
	if(r >1){
		new.beta.j <- beta.j.star
		accept <- accept+1
	}
	if(r <= 1){
		ifelse(u < r,new.beta.j <- beta.j.star,new.beta.j <- current.beta.j)
		ifelse(u < r,accept <- accept+1,accept <- 0)
	}
	
	out.beta.j <- list(new.beta.j=new.beta.j,accept.beta.j=accept)
	return(out.beta.j)
	
}



sample.mu <- function(m,beta.bar,Sigma,Lambda0,mu0,p){
	
	Lambda.n <- solve(solve(Lambda0)+(m*solve(Sigma)))
	mu.n <- Lambda.n%*%((solve(Lambda0)%*%matrix(mu0,p,1))+(solve(Sigma)%*%matrix(beta.bar,p,1)))
	
	new.mu <- mvrnorm(1,mu.n,Lambda.n)
	return(new.mu)
}



sample.Sigma <- function(beta,mu,m,eta0,S0,p){
	
	eta.n <- m+eta0
	S.mu <- matrix(0,p,p)
	for(j in 1:m){
		S.mu <- S.mu+(matrix(beta[j,]-mu,p,1)%*%t(matrix(beta[j,]-mu,p,1)))	
	}
	S.n <- S0+S.mu
	
	new.Sigma <- solve(rwish(eta.n,solve(S.n)))
	return(new.Sigma)
}


########################
### Initial values
########################

init.beta <- beta.ols
init.mu <- mu0
init.Sigma <- Lambda0



##############################
### Gibbs sampling algorithm
##############################


S <- 5000 # To save time, I changed from 50000 in lecture slides to 5000 here.

## We store the mouse specific regression coefficients in an array 
beta.MCMC <- array(0,dim=c(m,p,S))
## We store the mu and the elements of Sigma in a matrix
other.pars.MCMC <- matrix(0,S,ncol=30)
## This is for a new mouse yet to be sampled
beta.tilde.MCMC <- matrix(0,S,p)

## exp(mu*X) 
exp.mu.X.MCMC <- matrix(0,S,20)
## exp(betatilde*X) 
exp.betatilde.X.MCMC <- matrix(0,S,20)
## Ytilde 
Ytilde.MCMC <- matrix(0,S,20)

accept.beta <- rep(0,m)

set.seed(0) #The seed here is different from that in the lecture slides. Compare the results.
# let delta=0.25 for jumping distribution
for(k in 1:S){
	

	if(k==1){
		beta <- init.beta
		mu <- init.mu
		Sigma <- init.Sigma
	}
	beta.bar <- apply(beta,2,mean)

	new.mu <- sample.mu(m,beta.bar,Sigma,Lambda0,mu0,p)
	new.Sigma <- sample.Sigma(beta,new.mu,m,eta0,S0,p)
	
	new.beta <- matrix(0,m,p)
	for(j in 1:m){
		out.beta.j <- sample.beta.j(beta[j,],0.25,X,Y[j,],new.mu,new.Sigma,n.j)
		new.beta[j,] <- out.beta.j$new.beta.j
		accept.beta[j] <- accept.beta[j]+out.beta.j$accept.beta.j
	}
	
	
	mu <- new.mu
	Sigma <- new.Sigma
	beta <- new.beta
	
	beta.tilde.MCMC[k,] <- mvrnorm(1,new.mu,new.Sigma)
	exp.mu.X.MCMC[k,] <- exp(X%*%matrix(new.mu,p,1))
	exp.betatilde.X.MCMC[k,] <- exp(X%*%matrix(beta.tilde.MCMC[k,],p,1))
	Ytilde.MCMC[k,] <- rpois(20,exp.betatilde.X.MCMC[k,])
	
	
	beta.MCMC[,,k] <- beta
	other.pars.MCMC[k,1:5] <- mu
	other.pars.MCMC[k,6:30] <- as.numeric(Sigma)
	print(k)
	
}

### Acceptance rate
accept.beta/S


### Plots 

burnin <- 1000
thin <- 10


par(mfrow=c(1,3))

post.mean.mu.X <- apply(exp.mu.X.MCMC[seq((burnin+1),S,by=thin),],2,mean)
post.l95.mu.X <- apply(exp.mu.X.MCMC[seq((burnin+1),S,by=thin),],2,quantile,0.025)
post.u95.mu.X <- apply(exp.mu.X.MCMC[seq((burnin+1),S,by=thin),],2,quantile,0.975)

plot(seq(1:20),post.u95.mu.X,col="grey",lwd=2,lty=1,xlab="Location",ylab="Number of tumors",type="l",ylim=c(0,25))
lines(seq(1:20),post.l95.mu.X,col="grey",lwd=2,lty=1)
lines(seq(1:20),post.mean.mu.X,col="black",lwd=2,lty=1)


post.mean.betatilde.X <- apply(exp.betatilde.X.MCMC[seq((burnin+1),S,by=thin),],2,mean)
post.l95.betatilde.X <- apply(exp.betatilde.X.MCMC[seq((burnin+1),S,by=thin),],2,quantile,0.025)
post.u95.betatilde.X <- apply(exp.betatilde.X.MCMC[seq((burnin+1),S,by=thin),],2,quantile,0.975)

plot(seq(1:20),post.u95.betatilde.X,col="grey",lwd=2,lty=1,xlab="Location",ylab="Number of tumors",type="l",ylim=c(0,25))
lines(seq(1:20),post.l95.betatilde.X,col="grey",lwd=2,lty=1)
lines(seq(1:20),post.mean.betatilde.X,col="black",lwd=2,lty=1)



post.mean.Ytilde <- apply(Ytilde.MCMC[seq((burnin+1),S,by=thin),],2,mean)
post.l95.Ytilde <- apply(Ytilde.MCMC[seq((burnin+1),S,by=thin),],2,quantile,0.025)
post.u95.Ytilde <- apply(Ytilde.MCMC[seq((burnin+1),S,by=thin),],2,quantile,0.975)

plot(seq(1:20),post.u95.Ytilde,col="grey",lwd=2,lty=1,xlab="Location",ylab="Number of tumors",type="l",ylim=c(0,25))
lines(seq(1:20),post.l95.Ytilde,col="grey",lwd=2,lty=1)
lines(seq(1:20),post.mean.Ytilde,col="black",lwd=2,lty=1)



