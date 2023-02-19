#####  Metropolis algorithm for Poisson regression model
#####  Data: y_i| theta_i ~ Poisson(theta_i) with 
#####        E(Y_i | x_i)=theta_i=exp(beta_1 + beta_2 x_i + beta_3 x^2_i)
#####  Prior: beta=(beta_1,beta_2,beta_3) ~ N_3(0, 100*I_3)
#####  Proposal distribution: J_k(beta^star| beta^(k-1))=N_3(beta^(k-1), V)
#####                         where V=sigma2.hat*(X^prime*X)^{-1} and
#####                         sigma2.hat=variance of log(y_1+0.5),log(y_2+0.5),...,log(y_n+0.5)


library(MASS)
library(coda)
library(mvtnorm)

## Data

x.y.data <- read.table("../lecture_code/Sparrow.txt",header=TRUE)

x.y.data[1:3,]
y <- x.y.data$offspring
n <- length(y)
p <- 3
X <- matrix(cbind(x.y.data$intercept,x.y.data$age,x.y.data$age.squared),nrow=n,ncol=p)

## Prior 
beta0 <- rep(0,p)
Sigma <- 100*diag(1,3)

## Proposal variance
sigma2.hat <- var(log(y+0.5))
V <- sigma2.hat*solve((t(X)%*%X)) # cov mt of gaussian distribution


## Metropolis algorithm

sample.beta.metrop <- function(beta.k,var.proposal,y,X,beta0,Sigma){
	
	# this is a switch to keep track of the accepted values
	accept <- 0
	# proposing a candidate value
	beta.star <- mvrnorm(1,beta.k,var.proposal)
	
	# evaluating the metropolis ratio r=p(beta.star|y,X)/p(beta.k|y,X)
	# numerator: p(beta.star|y,X)=p(y_1,..,y_n|beta.star,X)*p(beta.star)/p(y|X), where p(y|X) is cancelled off.
	log.num <- sum(dpois(y,exp(X%*%beta.star),log=TRUE))+dmvnorm(beta.star,beta0,Sigma,log=TRUE)
	
	# denominator: p(beta.k|y,X)=p(y_1,..,y_n|beta.k,X)*p(beta.k)/p(y|X), where p(y|X) is cancelled off.
	log.den <- sum(dpois(y,exp(X%*%beta.k),log=TRUE))+dmvnorm(beta.k,beta0,Sigma,log=TRUE)
	
	r <- exp(log.num-log.den)
	
	# drawing a u~Unif(0,1)
	u <- runif(1,0,1)        # we can use bernoulli to replace this part
	
	# deciding whether to accept or reject
	if(r >1){
		new.beta <- beta.star	
	}
	if(r <=1){
		if(u <r){
			new.beta <- beta.star	
		}
		if(u >= r){
			new.beta <- beta.k	
		}
	}
	
	# Here we keep track of the number of acceptance 
	if(sum(new.beta==beta.star)==length(beta.star)){
		accept <- 1 
	}
	
	# I want to return two objects: the new value of beta and the acceptance counter.
	# We can return different type of outputs in a list format in R 
	out <- list(new.beta=new.beta,accept=accept)
	
	return(out)
}





### Running the Metropolis algorithm

S <- 10000


## Initial value for beta
glm.sparrow <- glm(y~x.y.data$age+x.y.data$age.squared,family=poisson)
beta.init <- as.numeric(glm.sparrow$coeff)


## We store the output in a matrix with S rows and 3 columns, 
beta.MCMC <- matrix(0,nrow=S,ncol=3)

no.accept <- 0

set.seed(0) # Note that the seed here is different from that for the lecture slides. Compare the results.

for(k in 1:S){
	
		if(k==1){
			beta <- beta.init	
		}
		out.beta <- sample.beta.metrop(beta,V,y,X,beta0,Sigma)
		new.beta <- out.beta$new.beta
	
		no.accept <- no.accept+out.beta$accept
		
		beta <- new.beta
		beta.MCMC[k,] <- beta
		print(k)
}


## Looking at the acceptance rate
accept.rate <- no.accept/S
accept.rate


## Trace plots
par(mfrow=c(1,3))
for(j in 1:3){
	plot(beta.MCMC[,j],ylab=paste("Beta_",j,sep=""),xlab="Iteration",main=paste("Trace plot for beta_",j,sep=""),type="l",lwd=1,lty=1)
}

# plot the first 500 iterations
par(mfrow=c(1,3))
for(j in 1:3){
  plot(beta.MCMC[1:500,j],ylab=paste("Beta_",j,sep=""),xlab="Iteration",main=paste("Trace plot for beta_",j,sep=""),type="l",lwd=1,lty=1)
}


# let burn-in period = 1000 iterations
burnin <- 1000 


## Autocorrelation function
par(mfrow=c(1,3))
for(j in 1:3){
	acf(beta.MCMC[(burnin+1):S,j],main=paste("ACF for beta_",j,sep=""),lag.max=40)
}


## Effective sample size
effectiveSize(beta.MCMC[(burnin+1):S,1])
effectiveSize(beta.MCMC[(burnin+1):S,2])
effectiveSize(beta.MCMC[(burnin+1):S,3])


## Thinning
thin.lag <- 20
sample.index <- seq((burnin+1),S,by=thin.lag)
length(sample.index)

## recheck the mixing of the sampled values
par(mfrow=c(1,3))
for(j in 1:3){
  acf(beta.MCMC[sample.index ,j],main=paste("ACF for beta_",j,sep=""),lag.max=40)
}

effectiveSize(beta.MCMC[sample.index,1])
effectiveSize(beta.MCMC[sample.index,2])
effectiveSize(beta.MCMC[sample.index,3])


## Marginal posterior densities
par(mfrow=c(1,3))
for(j in 1:3){
  hist(beta.MCMC[sample.index,j],xlab=paste("Beta_",j,sep=""),ylab="Density",
       main=paste("Marginal posterior density \n for beta_",j,sep=""),col="grey")
  abline(v=beta.init[j],col="red")
}

## Predicting theta_i

Xpred <- matrix(cbind(rep(1,6),seq(1:6),(seq(1:6)^2)),nrow=6,ncol=3)
theta.pred <- matrix(0,S,6)

for(k in 1:S){
	theta.pred[k,] <- exp(Xpred%*%beta.MCMC[k,])	
}


quantile.theta.pred <- apply(theta.pred[sample.index,],2,quantile,c(0.025,0.5,0.975))


par(mfrow=c(1,1))
plot(seq(1:6),quantile.theta.pred[2,],xlab="Age",ylab="Number of offspring",type="l",lwd=2,lty=1,col="black",ylim=range(as.numeric(theta.pred)),
main="Expected number of offspring as a function of age")
lines(seq(1:6),quantile.theta.pred[1,],col="red",lwd=2,lty=1)
lines(seq(1:6),quantile.theta.pred[3,],col="red",lwd=2,lty=1)




