#####  Metropolis algorithm for normal data with known variance
#####  Data: y_1,..,y_n | theta ~ N(theta, sigma^2) with sigma^2 known
#####  Prior: theta ~ N(mu, tau^2)
#####  Proposal distribution: J_k(theta^star|theta^(k-1))=N(theta^(k-1), delta^2)
#####  Posterior distribution: N(mu_n, tau^2_n)
#####  with:    tau^2_n = 1/((n/sigma^2)+(1/tau^2)) and mu_n=tau^2_n*((n/sigma^2)*ybar+(1/tau^2)*mu)



## Data

y <- c(9.37,10.18,9.16,11.60,10.33)
sigma2 <- 1
n <- length(y)
ybar <- mean(y)

## Prior 
mu <- 5
tau2 <- 10

## Posterior parameters
tau2.n <- 1/((n/sigma2)+(1/tau2))
mu.n <- tau2.n*(((n/sigma2)*ybar)+(mu/tau2))


## Metropolis algorithm

sample.theta.metrop <- function(theta.k,var.proposal,y,sigma2,mu,tau2){
	
	# this is a switch to keep track of the accepted values
	accept <- 0
	# proposing a candidate value
	theta.star <- rnorm(1,theta.k,sqrt(var.proposal))
	
	# evaluating the metropolis ratio r=p(theta.star|y)/p(theta.k|y)
	# numerator: p(theta.star|y)=p(y_1,..,y_n|theta.star)*p(theta.star)/p(y), where p(y) is cancelled off.
	log.num <- sum(dnorm(y,theta.star,sqrt(sigma2),log=TRUE))+dnorm(theta.star,mu,sqrt(tau2),log=TRUE)
	
	# denominator: p(theta.k|y)=p(y_1,..,y_n|theta.k)*p(theta.k)/p(y), where p(y) is cancelled off.
	log.den <- sum(dnorm(y,theta.k,sqrt(sigma2),log=TRUE))+dnorm(theta.k,mu,sqrt(tau2),log=TRUE)
	
	r <- exp(log.num-log.den)
	
	# drawing a u~Unif(0,1)
	u <- runif(1,0,1)
	
	# deciding whether to accept or reject
	if(r >1){
		new.theta <- theta.star	
	}
	if(r <=1){
		if(u <r){
			new.theta <- theta.star	
		}
		if(u >= r){
			new.theta <- theta.k	
		}
	}
	
	# Here we keep track of the number of acceptance 
	if(new.theta==theta.star){
		accept <- 1 
	}
	
	# I want to return two objects: the new value of theta and the acceptance counter.
	# We can return different type of outputs in a list format in R 
	out <- list(new.theta=new.theta,accept=accept)
	
	return(out)
}





### Running the Metropolis algorithm

S <- 10000


## Initial value for theta
theta.init <- 0

## Different values for the variance proposal delta^2: 1/32,1/2,2,32,64
delta2 <- c(1/32,1/2,2,32,64)

## We store the output in a matrix with S rows and 5 columns, one for each value of delta^2, the variance of
## the proposal distribution
theta.MCMC <- matrix(0,nrow=S,ncol=5)
## We store the acceptance rates in a 5x1 vector, one for each choice of delta^2
accept.rate <- rep(0,5)

## The outer loop is for the various values of delta^2; the inner loop are the MCMC iterations
set.seed(0) # Note that the seed here is different from that for the lecture slides. Compare the results.

for(j in 1:5){
	no.accept <- 0
	for(k in 1:S){
	
		if(k==1){
			theta <- theta.init	
		}
		out.theta <- sample.theta.metrop(theta,delta2[j],y,sigma2,mu,tau2)
		new.theta <- out.theta$new.theta
		no.accept <- no.accept+out.theta$accept
		
		theta <- new.theta
		theta.MCMC[k,j] <- theta
		print(k)
	}
	accept.rate[j] <- no.accept/S
}


## Looking at the acceptance rate
accept.rate


## Trace plots
par(mfrow=c(3,4))
for(j in 1:5){
	plot(theta.MCMC[,j],ylab="Theta",xlab="Iteration",main=paste("Trace plot for theta \n when delta^2 =",delta2[j]),type="l",lwd=1,lty=1)
	hist(theta.MCMC[,j],col="grey",prob=TRUE,xlab="Theta",ylab="density",main=paste("Histogram of MCMC samples for theta \n when delta^2 =",delta2[j]),breaks=100)
	lines(seq(min(theta.MCMC[,j]),max(theta.MCMC[,j]),length=1000),
		  dnorm(seq(min(theta.MCMC[,j]),max(theta.MCMC[,j]),length=1000),mu.n,sqrt(tau2.n)),col="red",type="l",lwd=2)
}


## Looking at the trace plots more closely
par(mfrow=c(1,5))
for(j in 1:5){
	plot(theta.MCMC[1:500,j],ylab="Theta",xlab="Iterations 1:500",main=paste("Trace plot for theta \n when delta^2 =",delta2[j]),type="l",lwd=1,lty=1)
}


## Autocorrelation function
# In practice, for the purpose of thinning, you should plot the acf plot after removing burn-in period.
par(mfrow=c(1,5))
for(j in 1:5){
	acf(theta.MCMC[,j],main=paste("ACF for theta \n when delta^2 =",delta2[j]),lag.max=100)
}


## Computing the lag-1 autocorrelation
lag1.acf <- rep(0,5)
for(j in 1:5){
	lag1.acf[j] <- acf(theta.MCMC[,j],lag.max=1,plot=FALSE)$acf[2,1,1]
}

lag1.acf






