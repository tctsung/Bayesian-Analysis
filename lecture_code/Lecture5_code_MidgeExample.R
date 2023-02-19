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



### trace plots
plot(phi[,1],xlab="Iteration number",ylab="Mu",main="Trace plot for Mu", type="l")
plot(phi[,2],xlab="Iteration number",ylab="Sigma^2",main="Trace plot for Sigma^2", type="l")

#only for the first 1000 iterations
plot(phi[1:1000,1],xlab="Iteration number",ylab="Mu",main="Trace plot for Mu", type="l")
plot(phi[1:1000,2],xlab="Iteration number",ylab="Sigma^2",main="Trace plot for Sigma^2", type="l")


###  acf plots
acf(phi[,1],lag=40,main="Autocorrelation function for Mu")
acf(phi[,2],lag=40,main="Autocorrelation function for Sigma^2")


burnin=500

### posterior distribution plots after 500 burnin iterations
plot(density(phi[(burnin+1):S,1]),xlab="Mu",main="Density estimate of the posterior \n marginal distribution of Mu", type="l")
plot(density(phi[(burnin+1):S,2]),xlab="Sigma^2",main="Density estimate of the posterior \n marginal distribution of Sigma^2", type="l")


### posterior summary after removing the 500 burnin iterations
## posterior means
mean(phi[(burnin+1):S,1])
mean(phi[(burnin+1):S,2])

## posterior medians
quantile(phi[(burnin+1):S,1],0.5)
quantile(phi[(burnin+1):S,2],0.5)

## 95% credible intervals based on posterior samples
quantile(phi[(burnin+1):S,1],c(0.025,0.975))
quantile(phi[(burnin+1):S,2],c(0.025,0.975))

## posterior SDs
sd(phi[(burnin+1):S,1])
sd(phi[(burnin+1):S,2])


