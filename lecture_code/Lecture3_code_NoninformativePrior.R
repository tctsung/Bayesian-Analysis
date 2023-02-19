#### Inference on a N(mu, sigma^2) when both mu and sigma^2 are unknown
#### using a non-informative improper prior 1/sigma^2
library(MASS) #for the function kde2d

y <- c(1.64, 1.70, 1.72, 1.74, 1.82, 1.82, 1.82, 1.90, 2.08)


ybar <- mean(y)
ybar
s2 <- var(y)
s2


### Sampling from the joint posterior distribution of mu and sigma2

B <- 10000

sigma2.samples <- rep(0,B)
mu.samples <- rep(0,B)

n <- length(y)

post.shape <- (n-1)/2
post.scale <- ((n-1)*s2)/2

set.seed(0) 
# set the random seed for reproducible work! 
# Note that I used a different seed for the lecture note. Compare the results.
for(i in 1:B){
	
	sigma2.samples[i] <- 1/rgamma(1,post.shape,post.scale)
	mu.samples[i] <- rnorm(1,ybar,sqrt(sigma2.samples[i]/n))
	print(i)
	
}

### Monte Carlo estimates of posterior mean, posterior variance and 95% credible interval for sigma^2
mean(sigma2.samples)
var(sigma2.samples)
quantile(sigma2.samples,0.025)
quantile(sigma2.samples,0.975)

## Exact values from theoretical distribution for sigma2
post.mean.sigma2 <- ((n-1)/(n-3))*s2
post.mean.sigma2

post.var.sigma2 <- 2*(((n-1)^2)/(((n-3)^2)*(n-4)))*(s2^2)
post.var.sigma2


### Monte Carlo estimates of posterior mean, posterior variance and 95% credible interval for mu
mean(mu.samples)
var(mu.samples)
quantile(mu.samples,0.025)
quantile(mu.samples,0.975)

## Exact values from theoretical distribution for mu
ybar+qt(0.025,n-1)*sqrt(s2/n)
ybar+qt(0.975,n-1)*sqrt(s2/n)


joint.dens <- kde2d(mu.samples,sigma2.samples)
par(mfrow=c(2,1))
contour(joint.dens,xlab="Mu",ylab="Sigma^2",xlim=c(1.7,1.9),ylim=c(0,0.04),drawlabels=FALSE,lwd=3,main="Monte Carlo estimate of the joint posterior density")
plot(mu.samples,sigma2.samples,type="p",col="black",pch=20,xlab="Mu",ylab="Sigma^2",main="Samples from the joint posterior distribution")



par(mfrow=c(1,1))
plot(density(sigma2.samples),xlab="Sigma^2",ylab="p(sigma^2 | y_1,..,y_n)",type="l",main="",lwd=3,lty=1,col="black")
abline(v=mean(sigma2.samples),col="red",lty=1,lwd=2)
abline(v=quantile(sigma2.samples,0.025),col="blue",lty=3,lwd=2)
abline(v=quantile(sigma2.samples,0.975),col="blue",lty=3,lwd=2)


par(mfrow=c(1,1))
plot(density(mu.samples),xlab="Mu",ylab="p(mu | y_1,.., y_n)",type="l",lwd=3,lty=1,col="black",main="")
abline(v=mean(mu.samples),col="red",lty=1,lwd=2)
abline(v=quantile(mu.samples,0.025),col="blue",lty=3,lwd=2)
abline(v=quantile(mu.samples,0.975),col="blue",lty=3,lwd=2)


par(mfrow=c(1,1))
plot(density(mu.samples),xlab="Mu",ylab="p(mu | y_1,.., y_n)",type="l",lwd=3,lty=1,col="black",main="")
abline(v=mean(mu.samples),col="red",lty=1,lwd=2)
abline(v=quantile(mu.samples,0.025),col="blue",lty=3,lwd=2)
abline(v=quantile(mu.samples,0.975),col="blue",lty=3,lwd=2)
abline(v=ybar+qt(0.025,n-1)*sqrt(s2/n),col="green",lty=3,lwd=2)
abline(v=ybar+qt(0.975,n-1)*sqrt(s2/n),col="green",lty=3,lwd=2)

