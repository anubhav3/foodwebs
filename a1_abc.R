rm(list=ls())
source("code/abc.functions.r")
source("code/Plot.matrix.r")

foodweb <- "Benguela Pelagic"

for (foodweb in c("Benguela Pelagic","Broadstone Stream","Broom","Capinteria","Caricaie Lakes","Coachella",
                  "EcoWEB41","EcoWEB60","Grasslands","Mill Stream","Sierra Lakes","Skipwith Pond",
                  "Small Reef","Tuesday Lake","Ythan","Grasslands"))
  {

## prior distribution
num <- 5000
a <- runif(num,-4,4)
ai <- runif(num,-2,2)
aj <- runif(num,-2,2)
r.b <- runif(num,-4,-1)
tol = 0.01 #tolerance to accept or reject parameter set
starting_parameter_values <- data.frame(a=a, ai=ai, aj=aj, r.b=r.b)

## select the web to analyse and load it
web.to.analyse <- foodweb
load(paste("data/food webs/", web.to.analyse, ".web.Rdata", sep=""))


## Plot the real food web matrix
name = paste(c(web.to.analyse,'.pdf'),collapse = '')
#pdf(name,width=7,height=5)
par(mfrow=c(3,2))


## fit the parameters
abc.RH.web <- abc.ratio(all.web.info, starting_parameter_values=starting_parameter_values, tol=tol)

## look at the parameter values
abc.RH.web$pars
## and the accuracy
abc.RH.web$acc

##Web connectivity
conn.real.web = format(round(sum(all.web.info$predation.matrix)/dim(all.web.info$predation.matrix)[1]^2, 2), nsmall = 2)
model.real.web =  format(round(sum(abc.RH.web$web)/dim(all.web.info$predation.matrix)[1]^2, 2), nsmall = 2)

#Proportion of link correct
prop.link.corr = format(round(sum(abc.RH.web$web==1 & all.web.info$predation.matrix==1)/sum(abc.RH.web$web), 2), nsmall = 2)

##Plot the real food webs
Plot.matrix(all.web.info$predation.matrix,title=paste("C",conn.real.web,sep="="))
## Plot the fitted food web matrix

Plot.matrix(abc.RH.web$web,title=paste("C=",model.real.web," Prop corr=",prop.link.corr,sep=""))

#Plotting the posterior distribution
hist((abc.RH.web$post)$ai, main= "ai",xlab="")
abline(v = mean((abc.RH.web$post)$ai),col="red")
hist((abc.RH.web$post)$aj, main= "aj",xlab="")
abline(v = mean((abc.RH.web$post)$aj),col="red")
hist((abc.RH.web$post)$a, main= "a",xlab="")
abline(v = mean((abc.RH.web$post)$a),col="red")
hist((abc.RH.web$post)$r.b, main= "r.b",xlab="")
abline(v = mean((abc.RH.web$post)$r.b),col="red")
mtext(web.to.analyse,  cex=1, line=40)
dev.copy2pdf(file = name,width=8,height=11)
print(paste(foodweb,"done!"))
}
