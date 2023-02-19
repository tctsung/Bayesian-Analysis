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

### Sampling from the joint posterior predictive distribution of Ytilde1 and Ytilde2

B <- 10000

theta1.post.samples <- rep(0,B)
theta2.post.samples <- rep(0,B)
ytilde1.samples <- rep(0,B)
ytilde2.samples <- rep(0,B)
dtilde.samples <- rep(0,B)


set.seed(0) 
# set the random seed for reproducible work! 
# Note that I used a different seed for the lecture note. Compare the results.
for(i in 1:B){

	theta1.post.samples[i] <- rgamma(1,atilde.theta1,btilde.theta1)
	theta2.post.samples[i] <- rgamma(1,atilde.theta2,btilde.theta2)
	ytilde1.samples[i] <- rpois(1,theta1.post.samples[i])
	ytilde2.samples[i] <- rpois(1,theta2.post.samples[i])
	dtilde.samples[i] <- ytilde1.samples[i]-ytilde2.samples[i]
	print(i)

}


## This is to compute the predictive distribution that Ytilde1 is greater than Ytilde2
mean(ytilde1.samples > ytilde2.samples)


#hist_dtilde
hist(dtilde.samples,xlab="College degree - No college degree",prob=TRUE,col="grey",
main="Posterior predictive distribution of \n difference in number of children",breaks=50)



### Sampling from the posterior predictive distribution of Ytilde1

B <- 10000

theta1.post.samples <- rep(0,B)
ytilde1.samples <- rep(0,B)

set.seed(0)
for(i in 1:B){
	
	theta1.post.samples[i] <- rgamma(1,atilde.theta1,btilde.theta1)
	ytilde1.samples[i] <- rpois(1,theta1.post.samples[i])
	print(i)
	
}


#hist_y1tilde
plot(seq(0,max(ytilde1.samples),by=1),as.numeric(table(ytilde1.samples))/B,type="n",xlab="Number of children",
main="Posterior predictive distribution of number of children \n Women without college degree",ylab="p(y | data)")
for(i in 0:max(ytilde1.samples)){
  segments(i,0,i,as.numeric(table(ytilde1.samples)/B)[i+1],col="black",lwd=6)
}




