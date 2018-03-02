rm(list=ls())

source("code/Plot.matrix.r")
source("code/ADBM.functions.r")

## make starting values
num <- 2
a <- seq(-4, 4, length=num)
ai <- seq(-2, 2, length=num)
aj <- seq(-2, 2, length=num)
r.b <- seq(-4, -1, length=num)
starting_parameter_values <- expand.grid(a=a, ai=ai, aj=aj, r.b=r.b)

## select the web to analyse and load it
web.to.analyse <- "Benguela Pelagic"
load(paste("data/food webs/", web.to.analyse, ".web.Rdata", sep=""))

## Plot the real food web matrix
Plot.matrix(all.web.info$predation.matrix)

## fit the parameters
NM.RH.web <- NM.ratio(all.web.info, starting_parameter_values=starting_parameter_values)

## look at the parameter values
NM.RH.web$pars
## and the accuracy
NM.RH.web$power

## Plot the fitted food web matrix
Plot.matrix(NM.RH.web$web)
