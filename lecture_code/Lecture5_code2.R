###  Sampling model: y_1,..,y_n | mu, sigma2 ~ N(mu, sigma2)
###  Prior: p(mu,sigma2)=p(mu)*p(sigma2)=N(mu0,tau2_0)*InverseGamma(nu0/2,nu0*sigma2_0/2)
###  Full conditionals:
###  - p(mu | y1,..,yn, sigma2) = N(mu_n,tau2_n)
###  - p(sigma2 | y1,.., yn, mu) = IGamma(nu_n/2, nu_n*sigma2_n(mu)/2)
###  where
###  tau2_n=1/((1/tau2_0)+(n/sigma2))
###  mu_n=tau2_n*((mu_0/tau2_0)+(n*ybar/sigma2))
###  nu_n = n+nu_0
###  sigma2_n(mu)=(1/nu_n)*(((n-1)*s2)+(n*(ybar-mu)^2)+(nu_0*sigma2_0))


###  Data

y <- c(1.64,1.70,1.72,1.74,1.82,1.82,1.82,1.90,2.08)
n <- length(y)
ybar <- mean(y)
s2 <- var(y)

### Prior

mu0 <- 1.9
tau2.0 <- 0.95^2
nu0 <- 1
sigma2.0 <- 0.01

### Functions to sample from the full conditionals
sample.sigma2 <- function(mu,nu0,sigma2.0,ybar,s2,n){

	nu.n <- n+nu0
	sigma2.n <- (1/nu.n)*(((n-1)*s2)+(n*(ybar-mu)^2)+(nu0*sigma2.0))
	
	post.shape <- nu.n/2
	post.scale <- (nu.n*sigma2.n)/2
	new.sigma2 <- 1/rgamma(1,post.shape,post.scale)
	return(new.sigma2)
}

sample.mu <- function(sigma2,mu0,tau2.0,ybar,n){
	
	tau2.n <- 1/((1/tau2.0)+(n/sigma2))
	mu.n <- tau2.n*((mu0/tau2.0)+(n*ybar/sigma2))
	
	new.mu <- rnorm(1,mu.n,sqrt(tau2.n))
	return(new.mu)
}



###
### Gibbs sampling algorithm

S <- 10000
phi <- matrix(0,S,2)

set.seed(1)
for(i in 1:S){

	if(i==1){
		# initial values
		phi[1,] <- c(ybar,s2)
	}
	
	if(i>1){
		phi[i,1] <- sample.mu(phi[(i-1),2],mu0,tau2.0,ybar,n)
		phi[i,2] <- sample.sigma2(phi[i,1],nu0,sigma2.0,ybar,s2,n)
	}
}



###

par(mfrow=c(1,1))
plot(phi[,1],phi[,2],xlab="Mu",ylab="Sigma^2",main="Samples from the joint posterior distribution \n Gibbs sampling algorithm",
type="p",pch=20,col="black")


par(mfrow=c(1,1))
image.plot(mu.grid,sigma2.grid,approx.joint.post.grid,xlab="Mu",ylab="Sigma^2",
main="p_D(mu, sigma^2 | y_1,..,y_n) and \n samples from joint posterior via Gibbs sampling",
col=heat.colors(100),xlim=c(1.6,2.0),ylim=c(0.001,0.025))

points(phi[,1],phi[,2],type="p",col="black")


plot(density(phi[,1]),xlab="Mu",main="Marginal posterior distribution p(mu | y_1,..,y_n) \n via Gibbs sampling",
type="l",lwd=3,lty=1,col="black")


plot(density(phi[,2]),xlab="Sigma^2",main="Marginal posterior distribution p(sigma^2 | y_1,..,y_n) \n via Gibbs sampling",
type="l",lwd=3,lty=1,col="black")

### Posterior summaries

## posterior means
mean(phi[,1])
mean(phi[,2])

## posterior medians
quantile(phi[,1],0.5)
quantile(phi[,2],0.5)



## 95% credible intervals based on posterior samples
quantile(phi[,1],c(0.025,0.975))
quantile(phi[,2],c(0.025,0.975))


plot(density(phi[,1]),xlab="Mu",main="Marginal posterior distribution p(mu | y_1,..,y_n) \n via Gibbs sampling",
type="l",lwd=3,lty=1,col="black")
abline(v=mean(phi[,1]),col="red",lwd=2)
abline(v=quantile(phi[,1],c(0.025,0.975)),col="grey",lwd=2)
abline(v=quantile(phi[,1],0.5),col="blue",lwd=2)


plot(density(phi[,2]),xlab="Sigma^2",main="Marginal posterior distribution p(sigma^2 | y_1,..,y_n) \n via Gibbs sampling",
type="l",lwd=3,lty=1,col="black")
abline(v=mean(phi[,2]),col="red",lwd=2)
abline(v=quantile(phi[,2],c(0.025,0.975)),col="grey",lwd=2)
abline(v=quantile(phi[,2],0.5),col="blue",lwd=2)


## posterior distribution of sigma

## posterior mean
mean(sqrt(phi[,2]))

## posterior median
quantile(sqrt(phi[,2]),0.5)

## 95% credible interval based on posterior samples
quantile(sqrt(phi[,2]),c(0.025,0.975))


plot(density(sqrt(phi[,2])),xlab="Sigma",main="Marginal posterior distribution p(sigma | y_1,..,y_n) \n via Gibbs sampling",
type="l",lwd=3,lty=1,col="black")
abline(v=mean(sqrt(phi[,2])),col="red",lwd=2)
abline(v=quantile(sqrt(phi[,2]),c(0.025,0.975)),col="grey",lwd=2)
abline(v=quantile(sqrt(phi[,2]),0.5),col="blue",lwd=2)


