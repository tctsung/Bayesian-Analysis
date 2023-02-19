library(coda)
library(MCMCpack)

tumor.data <- read.table("Tumor_counts.txt",header=FALSE)
Y <- as.matrix(tumor.data)

## Each row is a mouse
dim(tumor.data)

## Number of observations per mouse
n.j <- dim(tumor.data)[2]
## Number of mice
m <- dim(tumor.data)[1]


## Plotting the data
par(mfrow=c(1,1))

plot(1:20,Y[1,],type="l",xlab="Location",ylab="Number of tumors",lwd=2,lty=1,col="grey",ylim=c(0,18))
for(j in 2:m){
  lines(1:20,Y[j,],type="l",lwd=2,lty=1,col="grey")
}
lines(1:20,apply(Y,2,mean),col="black",lwd=2,lty=1,type="l")


## Creating the X matrix: each tumor number is a function of the tumor location
## and we model the relationship using a fourth degree polynomial
p <- 5
X <- matrix(0,nrow=n.j,ncol=5)
for(i in 1:n.j){
  X[i,] <- c(1,(i/20),(i/20)^2,(i/20)^3,(i/20)^4)
}


########################
### Prior distributions
########################

## Fitting a separate linear regression on the log average number of tumor per mouse
beta.ols <- matrix(0,nrow=m,ncol=5)
for(j in 1:m){
  log.data.mouse <- rep(0,20)
  for(i in 1:20){
    log.data.mouse[i] <- log(Y[j,i]+(i/20))	
  }
  beta.ols[j,] <- solve(t(X)%*%X)%*%(t(X)%*%matrix(log.data.mouse,nrow=20,ncol=1))	
}

mu0 <- apply(beta.ols,2,mean) # or colMeans(beta.ols)
Lambda0 <- cov(beta.ols)

eta0 <- 7
S0 <- Lambda0


init.beta <- beta.ols
init.mu <- mu0
init.Sigma <- Lambda0


Y.vect <- as.vector(t(Y)) # vectorize the matrix. as.vector() does vecterization of a matrix column by column, so here I transpose Y
Y[1,]-Y.vect[1:20] # check

X.rep <- rep(X[,2],21) #repeate X[,2] (i.e., i/20, i=1,...,20) for 21 times

id <- rep(1:21,each=20)

Data <- data.frame(Y=Y.vect, X1=X.rep, X2=X.rep^2, X3=X.rep^3, X4=X.rep^4, id=id)

model <- MCMChpoisson(fixed=Y~X1+X2+X3+X4, random=~X1+X2+X3+X4, group="id",
                      data=Data, burnin=1000, mcmc=10000, thin=1,verbose=1,
                      seed=0, 
                      beta.start=init.mu, Vb.start=init.Sigma, 
                      mubeta=mu0, Vbeta=Lambda0,
                      r=eta0, R=S0/eta0)

# Note that for Hierarchical Poisson Linear model,  
# the MCMCpack uses the blocked Gibbs sampler of Chib and Carlin (1999), not Gibbs with Metropolis algorithm
# and does not provide the tuning option for the acceptance rate in the blocked Gibbs sampler.
# But for ordinary Poisson Linear model, MCMCpoisson has the "tune" argument for the Metropolis algorithm.


attributes(model$mcmc) # burnin samples are discarded.
summary(model$mcmc)
mcmc.beta.X1 <- as.numeric(model$mcmc[,"beta.X1"]) # change to numeric, other the index uses the iteration #.
acf(mcmc.beta.X1)


