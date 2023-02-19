####  Basic hierarchical normal model:
####  Data: y_ij, j=1,..,m, i=1,..,n_j
####  1st stage model: y_ij | theta_j,sigma2 ~ N(theta_j,sigma2) for j=1,..,m; i=1,..,n_j
####  2nd stage model: theta_j | mu, tau2 ~ N(mu, tau2) for j=1,..,m
####  Priors:          sigma2 ~ InverseGamma(nu_0/2,nu_0*sigma2_0/2)
####                   mu ~ N(mu_0,tau2_0)
####                   tau2 ~ InverseGamma(eta_0/2, eta_0 * tau2_0/2)

####  Gibbs sampling algorithm will iteratively sample from all the full conditionals 
####  (for forms of full conditionals see Lecture 8)

library(coda)

#############
####  Data
#############

math.data <- read.table("schools_math_scores.txt",header=FALSE)

dim(math.data)
## Looking at the first 3 rows of data
math.data[1:3,]

## The first column in the data is the school indicator, the second is the student math scores
school <- math.data[,1]
y <- math.data[,2]

## The command table makes a table that reports the number of observation per school along with the school indicator
table(school)
## This is the vector with the number n_j of observations per school
n.j <- as.numeric(table(school))
## Number of schools
m <- length(n.j)

### Sample averages and sample variance by school
y.bar.school <- rep(0,m)
s2.school <- rep(0,m)

for(i in 1:m){
	y.bar.school[i] <- mean(y[which(school==i)])
    s2.school[i] <- var(y[which(school==i)])
}

## Average of the school-specific mean and variances
mean.y.bar.school <- mean(y.bar.school)
mean.s2.school <- mean(s2.school)
var.y.bar.school <- var(y.bar.school)


########################
### Prior distributions
########################
nu.0 <- 1
sigma2.0 <- 100
eta.0 <- 1
tau2.0 <- 100
mu.0 <- 50
gamma2.0 <- 25


########################
### Sampling functions
########################

sample.theta.j <- function(n.j,sigma2,mu,ybar.j,tau2){
	
	post.var <- 1/((n.j/sigma2)+(1/tau2))
	post.mean <- post.var*(((n.j/sigma2)*ybar.j)+(mu/tau2))
	
	new.theta.j <- rnorm(1,post.mean,sqrt(post.var))
	return(new.theta.j)
	
}



sample.sigma2 <- function(n.j,y,theta.j,nu.0,sigma2.0,m){
	
	# here we create a vector of length equal to the number of observations
	# where we have n_1 times theta_1, n_2 times theta_2,..., n_m times theta_m
	theta.j.expanded <- NULL
	for(i in 1:m){
		theta.j.expanded <- c(theta.j.expanded,rep(theta.j[i],n.j[i]))	
	}
	
	nu.n <- nu.0+sum(n.j)
	sigma2.n <- (1/nu.n)*((nu.0*sigma2.0)+sum((y-theta.j.expanded)^2))
	
	new.sigma2 <- 1/rgamma(1,(nu.n/2),((nu.n*sigma2.n)/2))
	return(new.sigma2)
}



sample.mu <- function(theta.j,m,tau2,mu.0,gamma2.0){
		
	post.var <- 1/((m/tau2)+(1/gamma2.0))
	theta.bar <- mean(theta.j)
	post.mean <- post.var*(((m/tau2)*theta.bar)+(mu.0/gamma2.0))
	
	new.mu <- rnorm(1,post.mean,sqrt(post.var))
	return(new.mu)
}



sample.tau2 <- function(theta.j,m,mu,eta.0,tau2.0){
	
	eta.n <- m+eta.0
	tau2.n <- (1/eta.n)*((eta.0*tau2.0)+sum((theta.j-mu)^2))
	
	new.tau2 <- 1/rgamma(1,(eta.n/2),((eta.n*tau2.n)/2))
	return(new.tau2)
}


########################
### Initial values
########################

init.theta.j <- y.bar.school
init.mu <- mean.y.bar.school
init.sigma2 <- mean.s2.school
init.tau2 <- var.y.bar.school



##############################
### Gibbs sampling algorithm
##############################


S <- 10000

## We store the school specific parameters, the theta_js in one matrix 
theta.MCMC <- matrix(0,nrow=S,ncol=m)
## We store the other parameters, mu, sigma2, tau2, in a second matrix
other.pars.MCMC <- matrix(0,S,ncol=3)
new.theta.j <- rep(0,m)

set.seed(0) #The random seed here is different from that for the slides of the lecture 8. Compare the results.
for(k in 1:S){
	

	if(k==1){
		theta.j <- init.theta.j
		mu <- init.mu
		sigma2 <- init.sigma2
		tau2 <- init.tau2
	}
	
	new.mu <- sample.mu(theta.j,m,tau2,mu.0,gamma2.0)
	new.tau2 <- sample.tau2(theta.j,m,new.mu,eta.0,tau2.0) 
	new.sigma2 <- sample.sigma2(n.j,y,theta.j,nu.0,sigma2.0,m,n)
		
	for(l in 1:m){
		new.theta.j[l] <- sample.theta.j(n.j[l],new.sigma2,new.mu,y.bar.school[l],new.tau2)	
	}
	mu <- new.mu
	tau2 <- new.tau2
	sigma2 <- new.sigma2
	theta.j <- new.theta.j
	
	theta.MCMC[k,] <- theta.j
	other.pars.MCMC[k,] <- c(mu,sigma2,tau2)
	print(k)
	
}



### MCMC diagnostic

## Trace plots
par(mfrow=c(2,2))
plot(other.pars.MCMC[,1],xlab="Iteration number",ylab="mu",type="l")
plot(other.pars.MCMC[,2],xlab="Iteration number",ylab="sigma2",type="l")
plot(other.pars.MCMC[,3],xlab="Iteration number",ylab="tau2",type="l")

par(mfrow=c(2,2))
plot(theta.MCMC[,1],xlab="Iteration number",ylab="theta_1",type="l")
plot(theta.MCMC[,2],xlab="Iteration number",ylab="theta_2",type="l")
plot(theta.MCMC[,3],xlab="Iteration number",ylab="theta_3",type="l")
plot(theta.MCMC[,4],xlab="Iteration number",ylab="theta_4",type="l")


## ACF plots
par(mfrow=c(2,2))
acf(other.pars.MCMC[,1],main="mu")
acf(other.pars.MCMC[,2],main="sigma2")
acf(other.pars.MCMC[,3],main="tau2")

par(mfrow=c(2,2))
acf(theta.MCMC[,1],main="theta_1")
acf(theta.MCMC[,2],main="theta_2")
acf(theta.MCMC[,3],main="theta_3")
acf(theta.MCMC[,4],main="theta_4")


## Effective sample size
effectiveSize(other.pars.MCMC[,1])
effectiveSize(other.pars.MCMC[,2])
effectiveSize(other.pars.MCMC[,3])

effectiveSize(theta.MCMC[,1])
effectiveSize(theta.MCMC[,2])
effectiveSize(theta.MCMC[,3])
effectiveSize(theta.MCMC[,4])

burnin <- 1000
L <- 1000
batch <- rep(seq(1:9),each=L)
## Box-plots of MCMC samples grouped in batches of L=1000 samples

par(mfrow=c(1,1))
boxplot(other.pars.MCMC[(burnin+1):S,1]~batch,ylab="Mu",xlab="Batch",col="grey")

par(mfrow=c(1,1))
boxplot(other.pars.MCMC[(burnin+1):S,2]~batch,ylab="Sigma2",xlab="Batch",col="grey")

par(mfrow=c(1,1))
boxplot(theta.MCMC[(burnin+1):S,3]~batch,ylab="Tau2",xlab="Batch",col="grey")

par(mfrow=c(1,1))
boxplot(theta.MCMC[(burnin+1):S,1]~batch,ylab="Theta_1",xlab="Batch",col="grey")

par(mfrow=c(1,1))
boxplot(theta.MCMC[(burnin+1):S,2]~batch,ylab="Theta_2",xlab="Batch",col="grey")


### Marginal posterior densities

plot(density(other.pars.MCMC[(burnin+1):S,1]),col="black",lwd=2,lty=1,xlab="Mu",ylab="Posterior density",main="Mu")
abline(v=quantile(other.pars.MCMC[(burnin+1):S,1],0.025),col="blue",lwd=2,lty=3)
abline(v=quantile(other.pars.MCMC[(burnin+1):S,1],0.975),col="blue",lwd=2,lty=3)
abline(v=median(other.pars.MCMC[(burnin+1):S,1]),col="purple",lwd=2,lty=1)
abline(v=mean(other.pars.MCMC[(burnin+1):S,1]),col="red",lwd=2,lty=1)


plot(density(other.pars.MCMC[(burnin+1):S,2]),col="black",lwd=2,lty=1,xlab="Sigma^2",ylab="Posterior density",main="Sigma^2")
abline(v=quantile(other.pars.MCMC[(burnin+1):S,2],0.025),col="blue",lwd=2,lty=3)
abline(v=quantile(other.pars.MCMC[(burnin+1):S,2],0.975),col="blue",lwd=2,lty=3)
abline(v=median(other.pars.MCMC[(burnin+1):S,2]),col="purple",lwd=2,lty=1)
abline(v=mean(other.pars.MCMC[(burnin+1):S,2]),col="red",lwd=2,lty=1)


plot(density(other.pars.MCMC[(burnin+1):S,3]),col="black",lwd=2,lty=1,xlab="Tau^2",ylab="Posterior density",main="Tau^2")
abline(v=quantile(other.pars.MCMC[(burnin+1):S,3],0.025),col="blue",lwd=2,lty=3)
abline(v=quantile(other.pars.MCMC[(burnin+1):S,3],0.975),col="blue",lwd=2,lty=3)
abline(v=median(other.pars.MCMC[(burnin+1):S,3]),col="purple",lwd=2,lty=1)
abline(v=mean(other.pars.MCMC[(burnin+1):S,3]),col="red",lwd=2,lty=1)


plot(density(theta.MCMC[(burnin+1):S,1]),col="black",lwd=2,lty=1,xlab="Theta_1",ylab="Posterior density",main="Theta_1")
abline(v=quantile(theta.MCMC[(burnin+1):S,1],0.025),col="blue",lwd=2,lty=3)
abline(v=quantile(theta.MCMC[(burnin+1):S,1],0.975),col="blue",lwd=2,lty=3)
abline(v=median(theta.MCMC[(burnin+1):S,1]),col="purple",lwd=2,lty=1)
abline(v=mean(theta.MCMC[(burnin+1):S,1]),col="red",lwd=2,lty=1)


plot(density(theta.MCMC[(burnin+1):S,2]),col="black",lwd=2,lty=1,xlab="Theta_2",ylab="Posterior density",main="Theta_2")
abline(v=quantile(theta.MCMC[(burnin+1):S,2],0.025),col="blue",lwd=2,lty=3)
abline(v=quantile(theta.MCMC[(burnin+1):S,2],0.975),col="blue",lwd=2,lty=3)
abline(v=median(theta.MCMC[(burnin+1):S,2]),col="purple",lwd=2,lty=1)
abline(v=mean(theta.MCMC[(burnin+1):S,2]),col="red",lwd=2,lty=1)


plot(density(theta.MCMC[(burnin+1):S,46]),col="black",lwd=2,lty=1,xlab="Theta_46",ylab="Posterior density",main="Theta_46")
abline(v=quantile(theta.MCMC[(burnin+1):S,46],0.025),col="blue",lwd=2,lty=3)
abline(v=quantile(theta.MCMC[(burnin+1):S,46],0.975),col="blue",lwd=2,lty=3)
abline(v=median(theta.MCMC[(burnin+1):S,46]),col="purple",lwd=2,lty=1)
abline(v=mean(theta.MCMC[(burnin+1):S,46]),col="red",lwd=2,lty=1)


plot(density(theta.MCMC[(burnin+1):S,82]),col="black",lwd=2,lty=1,xlab="Theta_82",ylab="Posterior density",main="Theta_82")
abline(v=quantile(theta.MCMC[(burnin+1):S,82],0.025),col="blue",lwd=2,lty=3)
abline(v=quantile(theta.MCMC[(burnin+1):S,82],0.975),col="blue",lwd=2,lty=3)
abline(v=median(theta.MCMC[(burnin+1):S,82]),col="purple",lwd=2,lty=1)
abline(v=mean(theta.MCMC[(burnin+1):S,82]),col="red",lwd=2,lty=1)


### Marginal posterior summaries

## Mu
mean(other.pars.MCMC[(burnin+1):S,1])
quantile(other.pars.MCMC[(burnin+1):S,1],c(0.025,0.975))
sd(other.pars.MCMC[(burnin+1):S,1])

## Sigma2
mean(other.pars.MCMC[(burnin+1):S,2])
quantile(other.pars.MCMC[(burnin+1):S,2],c(0.025,0.975))
sd(other.pars.MCMC[(burnin+1):S,2])


## Sigma
mean(sqrt(other.pars.MCMC[(burnin+1):S,2]))
quantile(sqrt(other.pars.MCMC[(burnin+1):S,2]),c(0.025,0.975))
sd(sqrt(other.pars.MCMC[(burnin+1):S,2]))

## Tau2
mean(other.pars.MCMC[(burnin+1):S,3])
quantile(other.pars.MCMC[(burnin+1):S,3],c(0.025,0.975))
sd(other.pars.MCMC[(burnin+1):S,3])


## Tau
mean(sqrt(other.pars.MCMC[(burnin+1):S,3]))
quantile(sqrt(other.pars.MCMC[(burnin+1):S,3]),c(0.025,0.975))
sd(sqrt(other.pars.MCMC[(burnin+1):S,3]))


## Theta
post.mean.theta <- apply(theta.MCMC[(burnin+1):S,],2,mean)



#### Shrinkage plots

plot(y.bar.school,post.mean.theta,xlab="School-specific sample average (ybar_j)",ylab="Posterior mean of theta_j",type="p",
col="grey",pch=20,main="Posterior mean of theta_j versus ybar_j")
abline(a=0,b=1,col="black")

plot(n.j,post.mean.theta-y.bar.school,xlab="Sample size for school j: n_j",ylab="Posterior mean of theta_j-ybar_j",type="p",
col="grey",pch=20,main="Posterior mean of theta_j minus ybar_j versus sample size n_j")
abline(h=0,col="black")


