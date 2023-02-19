

## 3D and 2D plots of the posterior probability that theta > some value w given the data in the 
## Beta-Binomial model

## The function that makes 3D plot is the persp function
## For help and documentation, type ?persp at the prompt in R
## The function persp takes as arguments: a x vector that gives the values on the x axis (this has
## to be in increasing order), a y vector with values on the y axis (also in increasing order),
## and a z matrix of dimension equal to the length of the x vector times the length of the y vector.
## The z matrix has at the (i,j) entry the value of the variable to be displayed on the third dimension
## corresponding to the i-th value of the x-vector and the j-th value of the y-vector, that is z[i,j] is 
## the value of the z variable at x[i] and y[j].

a <- seq(0,10,by=0.1)
b <- seq(0,10,by=0.1)
w <- 0.8
n <- 129
y <- 118


## Plotting the posterior distribution
#for(i in 1:length(a)){
#	for(j in 1:length(b)){
i=1
j=1
		plot(seq(0,1,by=0.01),dbeta(seq(0,1,by=0.01),a[i]+y,b[j]+(n-y)),type="l",col="black",lty=1,lwd=1,
			 xlim=c(0,1),xlab="Theta",ylab="P(theta | y)")
		title(main=paste("Posterior distribution a=",a[i]," b=",b[j],sep=""))
		abline(v=0.8,lty=3,col="red")
		
#	}
#}


post.prob.theta <- matrix(0,length(a),length(b))

for(i in 1:length(a)){
	for(j in 1:length(b)){
		
		post.a <- a[i]+y
		post.b <- b[j]+(n-y)
		
		post.prob.theta[i,j] <- 1-pbeta(w,post.a,post.b)
	}
}


persp(a,b,post.prob.theta,xlab="a",ylab="b",zlab="p(theta > w | y)")


# To rotate the figure, you can specify the value of the angles theta and phi by which to tilt the figure
persp(a,b,post.prob.theta,xlab="a",ylab="b",zlab="p(theta > 0.8 | y)",main="P(theta > 0.8 | y)",
theta=30,phi=20) # try ticktype = "detailed"




##### The 2D color plot may be better than the 3D plot.

## To make these 2D color plots it is necessary to download the R package fields
## If you do not know how to do this, please look at the reference R manual that is in th CTools website
## Once the package has been downloaded, you can install it by typing the command library(fields) at the prompt

library(fields)

## The function that makes 2D plot is the image.plot function
## For help and documentation, type ?image.plot at the prompt in R
## As the function persp, the function image.plot takes as arguments: a x vector that gives the values on the x axis (this has
## to be in increasing order), a y vector with values on the y axis (also in increasing order),
## and a z matrix of dimension equal to the length of the x vector times the length of the y vector.
## The z matrix has at the (i,j) entry the value of the variable to be displayed using color codes
## corresponding to the i-th value of the x-vector and the j-th value of the y-vector. 

a <- seq(0.1,10,by=0.1)
b <- seq(0.1,10,by=0.1)
w <- 0.8

n <- 129
y <- 118

post.prob.theta <- matrix(0,length(a),length(b))
post.mean.theta <- matrix(0,length(a),length(b))
theta0 <- matrix(0,length(a),length(b))
n0 <- matrix(0,length(a),length(b))

for(i in 1:length(a)){
	for(j in 1:length(b)){
		
		post.a <- a[i]+y
		post.b <- b[j]+(n-y)
		
		post.prob.theta[i,j] <- 1-pbeta(w,post.a,post.b)
		post.mean.theta[i,j] <- post.a/(post.a+post.b)
		theta0[i,j] <- a[i]/(a[i]+b[j])
		n0[i,j] <- a[i]+b[j]
	}
}



image.plot(a,b,post.prob.theta,xlab="a",ylab="b",main="p(theta > 0.8 | y)",col=heat.colors(100)[1:90])
contour(a,b,post.prob.theta, add=TRUE)

contour(a,b,post.prob.theta) # with no color
#https://stat.ethz.ch/R-manual/R-devel/library/graphics/html/contour.html

image.plot(n0,theta0,post.prob.theta,xlab="n0",ylab="theta0",main="p(theta > 0.8 | y)",col=heat.colors(100)[1:90])
#The contour() function does not work here, because n0 and theta0 are matrices and the 2D plot has some areas with no values. 

image.plot(n0,theta0,post.mean.theta,xlab="n0",ylab="theta0",main="E(theta |y)",col=heat.colors(100)[1:90])

