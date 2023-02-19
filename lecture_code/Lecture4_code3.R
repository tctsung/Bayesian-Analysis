### Model for the number of children for women without college degree and with college degree:
### the model is a Poisson-Gamma model with assumption of independence for the two priors

## prior for theta1
a.theta1 <- 2
b.theta1 <- 1


## prior for theta1
a.theta2 <- 2
b.theta2 <- 1


## data for women without college degree
n1 <- 111
sum.y1 <- 217

## data for women with college degree
n2 <- 44
sum.y2 <- 66

## posterior for theta1
atilde.theta1 <- a.theta1+sum.y1
btilde.theta1 <- b.theta1+n1

## posterior for theta2
atilde.theta2 <- a.theta2+sum.y2
btilde.theta2 <- b.theta2+n2



### Computing posterior p-values

B <- 10000
n1 <- 111
tobs <- 2

theta1.post.samples <- rep(0,B)
ytilde1.samples <- matrix(0,B,n1)
ttilde1.samples <- rep(0,B)

set.seed(0)
for(i in 1:B){
	
	theta1.post.samples[i] <- rgamma(1,atilde.theta1,btilde.theta1)
	ytilde1.samples[i,] <- rpois(n1,theta1.post.samples[i])
	ttilde1.samples[i] <- sum(ytilde1.samples[i,]==2)/sum(ytilde1.samples[i,]==1)
	print(i)
	
}

mean(ttilde1.samples > tobs)


#hist_ttilde1
hist(ttilde1.samples,xlab="Odds of having two children vs one",prob=TRUE,col="grey",
main="Posterior distribution of odds of having \n two children versus one",breaks=50)
abline(v=tobs)


