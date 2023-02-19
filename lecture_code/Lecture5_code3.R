####  MC approximation vs MCMC approximation to a distribution
####  Here we consider the approximation to the joint distribution p(delta, theta)
####  where delta is a discrete distribution that can take 3 values: 1, 2 or 3
####  p(delta=1)=0.45, p(delta=2)=0.1 and p(delta=3)=0.45
####  On the other hand, the conditional distribution of theta given delta is:
####  p(theta | delta=1)=N(theta; -3,1/3); p(theta| delta=2)=N(theta; 0, 1/3)
####  and p(theta | delta=3)=N(theta; 3, 1/3).
####  The marginal distribution of p(theta) is a mixture of the 3 normal distributions
####  with weights respectively equal to 0.45, 0.1 and 0.45



####  MC approximation to p(delta, theta)

B <- 1000

joint.samples.MC <- matrix(0,B,2)

set.seed(0)
for(i in 1:B){
	
	# sampling a value for delta from p(delta)
	joint.samples.MC[i,1] <- sample(c(1,2,3),1,replace=FALSE,prob=c(0.45,0.1,0.45))
	# sampling a value for theta from p(theta | delta) 
	if(joint.samples.MC[i,1]==1){
		joint.samples.MC[i,2] <- rnorm(1,-3,sqrt(1/3))
	}
	if(joint.samples.MC[i,1]==2){
		joint.samples.MC[i,2] <- rnorm(1,0,sqrt(1/3))
	}
	if(joint.samples.MC[i,1]==3){
		joint.samples.MC[i,2] <- rnorm(1,3,sqrt(1/3))
	}
}



hist(joint.samples.MC[,1],breaks=c(0.5,1.5,2.5,3.5),probability=TRUE,col="grey",xlab="Delta",main="MC approximation to p(delta)")

## this is to compute and plot the exact marginal distribution of theta
theta.values <- seq(-5,5,length=1000)

hist(joint.samples.MC[,2],col="grey",prob=TRUE,breaks=100,xlab="Theta",main="MC approximation to p(theta)")
lines(theta.values,0.45*dnorm(theta.values,-3,sqrt(1/3))+0.1*dnorm(theta.values,0,sqrt(1/3))+0.45*dnorm(theta.values,3,sqrt(1/3)),lty=1,lwd=2,col="red")




####  MCMC approximation to p(delta, theta): Gibbs sampling

## functions to sample from the full conditionals
sample.delta <- function(theta,mean.delta,sigma2.delta,prob.delta){

	cond.prob.delta <- rep(0,3)
	cond.prob.delta[1] <- prob.delta[1]*dnorm(theta,mean.delta[1],sqrt(sigma2.delta[1]))
	cond.prob.delta[2] <- prob.delta[2]*dnorm(theta,mean.delta[2],sqrt(sigma2.delta[2]))
	cond.prob.delta[3] <- prob.delta[3]*dnorm(theta,mean.delta[3],sqrt(sigma2.delta[3]))
	
	cond.prob.delta <- cond.prob.delta/sum(cond.prob.delta)
	
	new.delta <- sample(c(1,2,3),1,replace=FALSE,prob=cond.prob.delta)
	return(new.delta)
}


sample.theta <- function(delta,mean.delta,sigma2.delta){
	new.theta <- rnorm(1,mean.delta[delta],sqrt(sigma2.delta[delta]))
	return(new.theta)
}



### MCMC algorithm
S <- 1000

mean.delta <- c(-3,0,3)
sigma2.delta <- c(1/3,1/3,1/3)
prob.delta <- c(0.45,0.1,0.45)
joint.samples.MCMC <- matrix(0,S,2)

set.seed(1)
# using as initial values delta=1, theta=-3
init.delta <- sample(c(1,2,3),1,replace=FALSE,prob.delta)
init.theta <- rnorm(1,mean.delta[init.delta],1/3)

for(i in 1:S){
	
	if(i==1){
		joint.samples.MCMC[i,] <- c(init.delta,init.theta)
	}
	if(i>1){
		# sampling a value for delta from the full p(delta|theta)
		joint.samples.MCMC[i,1] <- sample.delta(joint.samples.MCMC[i-1,2],mean.delta,sigma2.delta,prob.delta)
		# sampling a value for theta from the full p(theta|delta)
		joint.samples.MCMC[i,2] <- sample.theta(joint.samples.MCMC[i,1],mean.delta,sigma2.delta)
	}
}



hist(joint.samples.MCMC[,1],breaks=c(0.5,1.5,2.5,3.5),probability=TRUE,col="grey",xlab="Delta",main="MCMC approximation to p(delta)")

hist(joint.samples.MCMC[,2],col="grey",prob=TRUE,breaks=100,xlab="Theta",main="MCMC approximation to p(theta)",xlim=c(-5,5))
lines(theta.values,0.45*dnorm(theta.values,-3,1/3)+0.1*dnorm(theta.values,0,1/3)+0.45*dnorm(theta.values,3,1/3),lty=1,lwd=2,col="red")


### Traceplot for delta and theta
plot(joint.samples.MCMC[,1],xlab="Iteration",main="Traceplot for delta",ylab="Delta")

plot(joint.samples.MCMC[,2],xlab="Iteration",main="Traceplot for theta",ylab="Theta")


##  Computing the autocorrelation
acf(joint.samples.MCMC[,2],plot=FALSE)

##  Computing the effective sample size
library(coda)
effectiveSize(joint.samples.MCMC[,2])

##  Determining the thinning parameter
acf(joint.samples.MCMC[,2],lag=150,main="Autocorrelation function of the MCMC samples for theta")


