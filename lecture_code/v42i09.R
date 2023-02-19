###################################################
### chunk number 1: 
###################################################
options(width=70, prompt="> ")


###################################################
### chunk number 2: killer-data
###################################################
library("MCMCpack") 
load("killamdt.rda") 
wilkerson <- 
MCMCregress(APRE1 ~ STRENGTH + COVER + WEAKMIN + NEWISSUE +
   SPONSOR + CHAMBER, data=killamdt, b0=0, B0=0.1, c0=2, d0=0.11,
   marginal.likelihood="Chib95") 


###################################################
### chunk number 3: killer-plot1
###################################################
plot(wilkerson[,1:2])


###################################################
### chunk number 4: Summary
###################################################
summary(wilkerson)  


###################################################
### chunk number 5: killer-plot2
###################################################
plot(density(wilkerson[,"(Intercept)"]+ wilkerson[,"SPONSOR"]),
xlim=c(0.2, 1.0), ylim=c(0,6.5),
        col="black", lwd=2, xlab="PRE", main="")
lines(density(wilkerson[,"(Intercept)"]), col="gray60", lwd=2)
legend(0.25, 6.6, legend=c("majority sponsored major weakening",
"minority sponsored major weakening"),
    col=c("black","gray60"), lwd=c(2,2), bty="n")


###################################################
### chunk number 6: models-killer
###################################################
model1 <- MCMCregress(APRE1 ~ STRENGTH + COVER, 
	data=killamdt, mcmc=10000, b0=0, B0=0.1, c0=2, d0=0.11, 
	marginal.likelihood="Chib95") 

model2 <- MCMCregress(APRE1 ~ STRENGTH + COVER + WEAKMIN +
    NEWISSUE, data=killamdt, mcmc=10000, b0=0, B0=0.1, c0=2,
    d0=0.11, marginal.likelihood="Chib95")

model3 <- MCMCregress(APRE1 ~ SPONSOR + CHAMBER, 
	data=killamdt, mcmc=10000, b0=0, B0=0.1, c0=2, d0=0.11, 
	marginal.likelihood="Chib95")

BF <- BayesFactor(model1, model2, model3, wilkerson)
summary(BF)


###################################################
### chunk number 7: mida1
###################################################

library("MCMCpack")
load("mida.rda")
plot(mida, type="h")


###################################################
### chunk number 8:  eval=FALSE
###################################################
## model1 <- MCMCpoissonChangepoint(mida, m=1, c0=13, d0=1,
##  	marginal.likelihood=c("Chib95"))    
## model2 <- MCMCpoissonChangepoint(mida, m=2, c0=13, d0=1,  
##  	marginal.likelihood=c("Chib95"))             
## model3 <- MCMCpoissonChangepoint(mida, m=3, c0=13, d0=1,  
##  	marginal.likelihood=c("Chib95")) 
## model4 <- MCMCpoissonChangepoint(mida, m=4, c0=13, d0=1,  
##  	marginal.likelihood=c("Chib95"))             
## model5 <- MCMCpoissonChangepoint(mida, m=5, c0=13, d0=1,  
##  	marginal.likelihood=c("Chib95"))                                          
## model6 <- MCMCpoissonChangepoint(mida, m=6, c0=13, d0=1,  
##  	marginal.likelihood=c("Chib95"))                                          
## model7 <- MCMCpoissonChangepoint(mida, m=7, c0=13, d0=1,  
##  	marginal.likelihood=c("Chib95"))                                          
## model8 <- MCMCpoissonChangepoint(mida, m=8, c0=13, d0=1,  
## 	marginal.likelihood=c("Chib95"))


###################################################
### chunk number 9: mida-model
###################################################
load("pcModels.rda")
BF <- BayesFactor(model1, model2, model3, model4, model5, model6, model7,
   model8)
summary(BF)


###################################################
### chunk number 10: mida2
###################################################

attr(model6, "y") <- ts(attr(model6, "y"), start = 1816)
plotState(model6, legend.control=c(1813, 0.75))  


###################################################
### chunk number 11: mida-fig3
###################################################

plotChangepoint(model6, main = "", verbose=FALSE)


