rm(list=ls())

source("code/Plot.matrix.r")
source("code/ADBM.functions.r")

web.to.analyse <- "Benguela Pelagic"

load(paste("data/food webs/", web.to.analyse, ".web.Rdata", sep=""))

Plot.matrix(all.web.info$predation.matrix)

NM.RH.web <- NM.ratio(all.web.info, init.par.file="data/parameterisation/ratio.initial.pars.txt")

Plot.matrix(NM.RH.web$web)

