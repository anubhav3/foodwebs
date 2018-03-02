rm(list=ls())

source("code/Plot.matrix.r")
source("code/ADBM.functions.r")

## make starting values
num <- 2
a <- 10^seq(-4, 4, length=num)
ai <- seq(-2, 2, length=num)
aj <- seq(-2, 2, length=num)
r.b <- 10^seq(-4, -1, length=num)
starting_parameter_values <- expand.grid(a=a, ai=ai, aj=aj, r.b=r.b)

web.to.analyse <- "Benguela Pelagic"

load(paste("data/food webs/", web.to.analyse, ".web.Rdata", sep=""))

Plot.matrix(all.web.info$predation.matrix)

NM.RH.web <- NM.ratio(all.web.info, starting_parameter_values=starting_parameter_values)

Plot.matrix(NM.RH.web_match1s$web)
Plot.matrix(NM.RH.web_accuracy$web)

sum(NM.RH.web_accuracy$web)
sum(NM.RH.web_match1s$web)
##NM.RH.web_match1s <- NM.RH.web

##NM.RH.web_accuracy <- NM.RH.web

NM.RH.web_match1s$pars
NM.RH.web_accuracy$pars
