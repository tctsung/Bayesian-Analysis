###  Sampling model: y_1,..,y_n | mu, sigma2 ~ N(mu, sigma2)
###  Prior: p(mu,sigma2)=p(mu)*p(sigma2)=N(mu0,tau2_0)*InverseGamma(nu0/2,nu0*sigma2_0/2)

###  Data

y <- c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
ybar <- mean(y)
s2 <- var(y)

### Prior

mu0 <- 1.9
tau2.0 <- 0.95^2
nu0 <- 1
sigma2.0 <- 0.01



####  Approximation of the joint posterior distribution

G <- 100
H <- 100

mu.grid <- seq(1.505,2.00,length=G)
sigma2.grid <- seq(0.001,0.1,length=H)


approx.joint.post.grid <- matrix(0,G,H)

### This is the function that evaluates the inverse gamma density (on the log scale)
log.inverse.gamma.dens <- function(x,a,b){
	out <- a*log(b)-log(gamma(a))-(a+1)*log(x)-(b/x)
	return(out)
}


### This is the function that approximates the joint posterior density on a grid via p_D(mu, sigma2 | y_1,..,y_n) 
### To avoid working with too small values, we evaluate all the probabilities on the log scale and then we
### transform back on the original scale at the end

for(j in 1:G){
	for(k in 1:H){
	
		# Likelihood on the log scale for mu equal to the j-th value on the grid for my and with sigma2 equal to the k-th value of the grid for sigma2 
		log.lik <- sum(dnorm(y,mu.grid[j],sqrt(sigma2.grid[k]),log=TRUE))
	
		# Prior for mu evaluated at the j-th value of the grid for mu on the log-scale
		log.prior.mu.grid <- dnorm(mu.grid[j],mu0,sqrt(tau2.0),log=TRUE)
		
		# Prior for sigma2 evaluated at the k-th value of the grid for sigma2 on the log-scale
		a <- nu0/2
		b <- nu0*sigma2.0/2
		log.prior.sigma2.grid <- log.inverse.gamma.dens(sigma2.grid[k],a,b)
		
		
		# This is the un-normalized approximate discrete joint posterior distribution
		approx.joint.post.grid[j,k] <- exp(log.lik+log.prior.mu.grid+log.prior.sigma2.grid)
	}
}


approx.joint.post.grid <- approx.joint.post.grid/sum(approx.joint.post.grid)


### This is to plot the approximate discrete joint posterior distribution
library(fields)

par(mfrow=c(1,1))
image.plot(mu.grid,sigma2.grid,approx.joint.post.grid,xlab="Mu",ylab="Sigma^2",main="p_D(mu, sigma^2 | y_1,..,y_n)",
col=heat.colors(100),xlim=c(1.6,2.0),ylim=c(0.001,0.025))



### To obtain the approximate discrete marginal posterior distribution of mu given y_1, y_2, .., y_n
approx.marg.post.mu.grid <- rep(0,G)
for(j in 1:G){
	approx.marg.post.mu.grid[j] <- sum(approx.joint.post.grid[j,])	
}

## Alternatively, we can obtain approx.marg.post.mu.grid simply adding up the elements of the matrix approx.joint.post.grid row by row
approx.marg.post.mu.grid.v2 <- apply(approx.joint.post.grid,1,sum)


### To obtain the approximate discrete marginal posterior distribution of sigma2 given y_1, y_2, .., y_n
approx.marg.post.sigma2.grid <- rep(0,H)
for(k in 1:H){
	approx.marg.post.sigma2.grid[k] <- sum(approx.joint.post.grid[,k])	
}

## Alternatively, we can obtain approx.marg.post.sigma2.grid simply adding up the elements of the matrix approx.joint.post.grid column by column
approx.marg.post.sigma2.grid.v2 <- apply(approx.joint.post.grid,2,sum)


### This is to plot the two approximate discrete marginal posterior distributions
par(mfrow=c(1,2),mai=c(0.8,0.8,0.5,0.5))
plot(mu.grid,approx.marg.post.mu.grid,xlab="Mu",ylab="p(mu | y_1,.., y_n)",type="l",lwd=2,lty=1,col="black")
plot(sigma2.grid,approx.marg.post.sigma2.grid,xlab="Sigma^2",ylab="p(sigma^2 | y_1,.., y_n)",type="l",lwd=2,lty=1,col="black")




#### Rejection sampling

### Example: here the density from which we want to sample is L(theta)*p(theta)=Beta(3,10)
### This density takes as maximum value 3.585
### Let M be 4 and g(theta) equal to the uniform density on [0,1]

theta.values <- seq(0,1,length=1000)
rawDensity <- dbeta(theta.values, 3,10)
plot(theta.values,rawDensity,type="l",lwd=1,lty=1,col="black",xlab="Theta",ylab="density",ylim=c(0,5))
lines(theta.values,rep(4,1000),type="l",lwd=1,lty=1,col="black")
text(0.58,1.3,"L(theta)*p(theta)",cex=1.3,col="red")
text(0.4,4.2,"M*g(theta)",cex=1.3,col="red")



### Rejection sampling algorithm
RejectionSampling <- function(M,n){

	accepted.samples <- NULL
	for(i in 1:n){
		ok <- 0 
		while(ok<1){
			# sampling a value from the envelope density
			theta.i <- runif(1,min=0,max=1)
			
			# sampling a uniform value from the uniform distribution
			u <- runif(1,min=0,max=1)
			
			# determining whether to accept or reject the sampled theta.i
			# if accepted, add to the sample of accepted thetas and change ok to 1
			if(M*u*1 <= dbeta(theta.i,3,10)){
				ok <- 1 
				accepted.samples <- c(accepted.samples,theta.i)
			}
			# if not accepted, draw another sample for theta.i from the envelope density and another value of u
		}
	}
	return(accepted.samples)
}


## Here we compare the histogram of the sample generated via rejection sampling
## with the density of a beta(3,10) random variable
set.seed(0)
simulatedDensity <- RejectionSampling(4,3000)

par(mfrow=c(1,1))
hist(simulatedDensity,xlab="theta",main="",breaks=50,col="grey",prob=TRUE)
lines(theta.values,rawDensity,col="red",lty=1,lwd=2)
