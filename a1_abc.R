rm(list=ls())

source("code/abc.functions.r")
source("code/Plot.matrix.r")

## prior distribution
num <- 5000
a <- runif(num,-4,4)
ai <- runif(num,-2,2)
aj <- runif(num,-2,2)
r.b <- runif(num,-4,-1)
tol = 0.05 #Percentage of samples to be chosen
starting_parameter_values <- data.frame(a=a, ai=ai, aj=aj, r.b=r.b, tol=tol)

## select the web to analyse and load it
web.to.analyse <- "Benguela Pelagic"
load(paste("data/food webs/", web.to.analyse, ".web.Rdata", sep=""))


## Plot the real food web matrix
name = paste(c(web.to.analyse,'.pdf'),collapse = '')
pdf(name,width=7,height=5)
par(mfrow=c(3,2))
Plot.matrix(all.web.info$predation.matrix)

## fit the parameters
abc.RH.web <- abc.ratio(all.web.info, starting_parameter_values=starting_parameter_values)

## look at the parameter values
abc.RH.web$pars
## and the accuracy
abc.RH.web$acc

## Plot the fitted food web matrix
Plot.matrix(abc.RH.web$web)

#Plotting the posterior distribution
hist((abc.RH.web$post)$ai, main= "ai",xlab="")
abline(v = mean((abc.RH.web$post)$ai),col="red")
hist((abc.RH.web$post)$aj, main= "aj",xlab="")
abline(v = mean((abc.RH.web$post)$aj),col="red")
hist((abc.RH.web$post)$a, main= "a",xlab="")
abline(v = mean((abc.RH.web$post)$a),col="red")
hist((abc.RH.web$post)$r.b, main= "r.b",xlab="")
abline(v = mean((abc.RH.web$post)$r.b),col="red")
mtext(web.to.analyse,  cex=1, line=26)
dev.off()