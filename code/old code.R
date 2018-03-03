
####################################################
## FOOD WEB MODEL FUNCTIONS
####################################################

Allometric.EHL <- function(M,
                           e,
                           h, hi, hj,
                           a, ai, aj,
                           n, ni=-3/4,
                           h.method=1){
  ## HANDLING TIMES
  if(h.method==1){
    get.h <- function(Mi, Mj, h, hi, hj)
      h * Mi^hi * Mj^hj
  }
  if(h.method==2){
    get.h <- function(Mi, Mj,
                      h1, h1i, h1j,
                      h2, h2i, h2j)
      h1 * Mi^h1i * Mj^h1j + h2 * Mi^h2i * Mj^h2j
  }
  if(h.method==3){
    get.h <- function(Mi, Mj, r.a, r.b)
      ifelse((r.b-Mi/Mj)>0,
             r.a/(r.b-Mi/Mj),
             Inf)
  }
  ## in matrix H resources are rows and consumers are columns
  H <- outer(M, M, get.h, h=h, hi=hi, hj=hj)
  
  ## ENCOUNTER RATES: consumer - resource mass specific encounter rates
  get.a <- function(Mi, Mj,
                    a, ai, aj)
    a * Mi^ai * Mj^aj
  A <- outer(M, M, get.a,
             a=a, ai=ai, aj=aj)
  L <- A* n*M^ni
  
  if(sum(order(M)==1:length(M))!=length(M))
    stop("Body sizes not sorted")
  
  ## energy values
  E <- e*M
  
  return = list(E=E, H=H, L=L)
  
}

Allometric.EHL.diet <- function(Mi, Mj,
                                e,
                                h, hi, hj,
                                a, ai, aj,
                                n, ni=-3/4,
                                h.method=1){
  
  if(sum(order(Mi)==1:length(Mi))!=length(Mi))
    stop("Body sizes not sorted")
  
  H <- h*Mi^hi*Mj^hj
  A = a * Mi^ai * Mj^aj
  L <- A* n*Mi^ni
  E <- e*Mi    
  list(E=E, H=H, L=L)
}

Ratio.allometric.EHL.diet <- function(Mi, Mj,
                                      e,
                                      r.a, r.b,
                                      a, ai, aj,
                                      n, ni=-3/4){
  
  if(sum(order(Mi)==1:length(Mi))!=length(Mi))
    stop("Body sizes not sorted")
  
  if(!r.b==0)
    H = ifelse((r.b-Mi/Mj)>0,
               r.a/(r.b-Mi/Mj),
               Inf)
  if(r.b==0)
    H = rep(r.a, length(Mi))
  
  
  A = a * Mi^ai * Mj^aj
  L <- A* n*Mi^ni
  
  ## energy values
  E <- e*Mi
  
  list(E=E, H=H, L=L)
  
}




Get.diet <- function(EHL){
  P <- EHL[[1]]/EHL[[2]]
  p <- P
  order.by.p <- order(p, decreasing=T)
  p <- p[order.by.p]
  Lj <- EHL[[3]][order.by.p]
  hj <- EHL[[2]][order.by.p]
  Ej <- EHL[[1]][order.by.p]
  cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum( Lj * hj))
  n.d <- which(cumulative.profit==max(cumulative.profit))
  diet = c(rep(T, n.d), rep(F, length(p) - n.d))[order(order.by.p)]
  diet
}

Get.diet.E <- function(EHL){
  P <- EHL[[1]]/EHL[[2]]
  p <- P
  order.by.p <- order(p, decreasing=T)
  p <- p[order.by.p]
  Lj <- EHL[[3]][order.by.p]
  hj <- EHL[[2]][order.by.p]
  Ej <- EHL[[1]][order.by.p]
  cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum( Lj * hj))
  max(cumulative.profit)
}


Get.web.old <- function(EHL, energy.intake=F){
  
  ##E <- EHL[[1]]
  ##H <- EHL[[2]]
  ##L <- EHL[[3]]
  S <- length(EHL[[1]])
  
  web <- matrix(0, S, S)
  overall.energy <- numeric(S)
  per.species.energy <- matrix(0, S, S)
  
  ## in matrix P, columns are cosumers and contain profit of that consumer
  ## feeding on each prey (row)
  P <- EHL[[1]]/EHL[[2]]
  
  ## split code depending on whether encounter rates are predator specific or not
  if(is.matrix(EHL[[3]])){
    for(j in 1:S){
      
      ## ordering of p required
      p <- P[,j]
      order.by.p <- order(p, decreasing=T)
      p <- p[order.by.p]
      Lj <- EHL[[3]][,j][order.by.p]
      hj <- EHL[[2]][,j][order.by.p]
      Ej <- EHL[[1]][order.by.p]
      
      cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum( Lj * hj))
      ##cumulative.profit[length(E)] <- NA
      
      if(!all(p==0)){
        web[,j] <- c(1, cumulative.profit[1:(length(EHL[[1]])-1)] < p[2:length(EHL[[1]])])[order(order.by.p)]
        overall.energy[j] <- cumulative.profit[sum(web[,j])]
      }
      energies <- c(Ej * Lj)[1:sum(web[,j])] / (1 + cumsum( Lj * hj)[sum(web[,j])]) 
      all.energies <- c(energies, rep(0, S-length(energies)))
      
      per.species.energy[,j] <- all.energies[order(order.by.p)]
      
    }
  }
  if(is.vector(EHL[[3]])){
    for(j in 1:S){
      
      ## ordering of p required
      p <- P[,j]
      order.by.p <- order(p, decreasing=T)
      p <- p[order.by.p]
      Lj <- EHL[[3]][order.by.p]
      hj <- EHL[[2]][,j][order.by.p]
      Ej <- EHL[[1]][order.by.p]
      
      cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum( Lj * hj))
      
      dj <- max(which(cumulative.profit==max(cumulative.profit)))
      web[,j] <- c(rep(1, dj), rep(0, S-dj))[order(order.by.p)]
      
      overall.energy[j] <- cumulative.profit[sum(web[,j])]    
      
      energies <- c(Ej * Lj)[1:sum(web[,j])] / (1 + cumsum( Lj * hj)[sum(web[,j])]) 
      all.energies <- c(energies, rep(0, S-length(energies)))
      
      per.species.energy[,j] <- all.energies[order(order.by.p)]
      
    }
  }
  if(energy.intake)
    result <- list(web=web, overall.flux=overall.energy, per.species.flux=per.species.energy)
  else
    result <- web
  result
}

Empirical.allometric.DBM <- function(all.web.info, parameters, handling.function, NITS=10, optimise.C=F){
  
  S <- length(all.web.info$species.sizes)
  M <- all.web.info$species.sizes
  target.C <- sum(all.web.info$predation.matrix) / S^2
  real.web <- all.web.info$predation.matrix
  
  if(parameters=="EAT" & all.web.info[[3]]=="length"){
    
    ## energy
    e = 1
    
    ## density
    if(.Platform$OS.type=="windows")
      load(paste(path, "data//type II size investigation//n.length.coefs.Rdata", sep=""))
    if(.Platform$OS.type=="unix")
      load("~/allometric.web/data/Type II size investigation//n.length.coefs.Rdata")
    
    m.n = n.length.coefs[1,1]
    m.ni = n.length.coefs[2,1]
    sd.n = n.length.coefs[1,2]
    sd.ni = n.length.coefs[2,2]
    
    ## handling times
    if(.Platform$OS.type=="unix")
      load("~/allometric.web/data/Type II size investigation//H.length.coefs.Rdata")
    m.h = H.length.coefs[1,1]
    m.hi = H.length.coefs[2,1]
    m.hj = H.length.coefs[3,1]
    sd.h = H.length.coefs[1,2]
    sd.hi = H.length.coefs[2,2]
    sd.hj = H.length.coefs[3,2]
    
    ## encounter rate
    if(.Platform$OS.type=="unix")
      load("~/allometric.web/data/Type II size investigation//a.length.coefs.Rdata")
    m.a = a.length.coefs[1,1]
    m.ai = a.length.coefs[2,1]
    m.aj = a.length.coefs[3,1]
    sd.a = a.length.coefs[1,2]
    sd.ai = a.length.coefs[2,2]
    sd.aj = a.length.coefs[3,2]
    
  }
  
  
  if(parameters=="EAT" & all.web.info[[3]]=="mass"){
    
    ## energy
    e = 1
    
    ## density
    if(.Platform$OS.type=="unix")
      load("~/allometric.web/data/Type II size investigation//n.mass.coefs.Rdata")
    m.n = 10^n.mass.coefs[1,1]
    m.ni = n.mass.coefs[2,1]
    sd.n = n.mass.coefs[1,2]
    sd.ni = n.mass.coefs[2,2]
    
    ## handling times
    if(.Platform$OS.type=="unix")
      load("~/allometric.web/data/Type II size investigation//H.mass.coefs.Rdata")
    m.h = 10^H.mass.coefs[1,1]
    m.hi = H.mass.coefs[2,1]
    m.hj = H.mass.coefs[3,1]
    sd.h = H.mass.coefs[1,2]
    sd.hi = H.mass.coefs[2,2]
    sd.hj = H.mass.coefs[3,2]
    
    ## encounter rate
    if(.Platform$OS.type=="unix")
      load("~/allometric.web/data/Type II size investigation//a.mass.coefs.Rdata")
    m.a = 10^a.mass.coefs[1,1]
    m.ai = a.mass.coefs[2,1]
    m.aj = a.mass.coefs[3,1]
    sd.a = a.mass.coefs[1,2]
    sd.ai = a.mass.coefs[2,2]
    sd.aj = a.mass.coefs[3,2]
    
  }
  
  if(parameters=="EWT" & all.web.info[[3]]=="mass"){        
    
    ## energy
    e = 1
    
    ## density
    ## density (per m2) from Duarte 1987 Oecologia
    #temp.x <- 10^c(0, 8.98)*1e-12
    #temp.y <- 10^c(8.53, 0)*1e6
    #lm(log10(temp.y) ~ log10(temp.x))
    m.n = 3.13  ## I've checked that this does not need to be logged
    sd.n = 0
    m.ni = -0.95
    sd.ni= 0
    
    if(.Platform$OS.type=="windows")
      T.p <- read.csv(paste(path, "data\\Thompson\\pars.mass.csv", sep="\\"), header=T)[,-1]
    if(.Platform$OS.type=="unix")
      T.p <- read.csv("~/allometric.web/data/Thompson/pars.mass.csv", header=T)[,-1]
    
    ## handling times
    m.h = T.p[1,1]
    m.hi = T.p[3,1]
    m.hj = T.p[2,1]
    sd.h = T.p[1,2]
    sd.hi = T.p[3,2]
    sd.hj = T.p[2,2]
    
    ## encounter rate
    m.a = T.p[4,1]
    m.ai = T.p[6,1]
    m.aj = T.p[5,1]
    sd.a = T.p[4,2]
    sd.ai = T.p[6,2]
    sd.aj = T.p[5,2]
  }
  
  if(parameters=="EWT" & all.web.info[[3]]=="length"){        
    
    ## energy
    e = 1
    
    ## density
    ## density (per m2) from Duarte 1987 Oecologia
    #temp.x <- 10^c(0, 8.98)*1e-12
    #temp.y <- 10^c(8.53, 0)*1e6
    #lm(log10(temp.y) ~ log10(temp.x))
    m.n = 3.13  ## I've checked that this does not need to be logged
    sd.n = 0
    m.ni = -0.95
    sd.ni = 0
    
    if(.Platform$OS.type=="windows")
      T.p <- read.csv(paste(path, "data\\Thompson\\pars.length.csv", sep="\\"), header=T)[,-1]
    if(.Platform$OS.type=="unix")
      T.p <- read.csv("~/allometric.web/data/Thompson/pars.length.csv", header=T)[,-1]
    
    ## handling times
    m.h = T.p[1,1]
    m.hi = T.p[3,1]
    m.hj = T.p[2,1]
    sd.h = T.p[1,2]
    sd.hi = T.p[3,2]
    sd.hj = T.p[2,2]
    
    ## encounter rate
    m.a = T.p[4,1]
    m.ai = T.p[6,1]
    m.aj = T.p[5,1]
    sd.a = T.p[4,2]
    sd.ai = T.p[6,2]
    sd.aj = T.p[5,2]
  }
  
  
  #performance <- list(1)
  #observed.connectance = numeric(length=NITS)
  #predicted.connectance =  numeric(length=NITS)
  #web.id <- numeric(length=NITS)
  #prop.correct <- numeric(length=NITS)
  webs <- list(1)
  EHLs <- list(1)
  best <- list(1)
  
  ## power parameters    
  M <- all.web.info$species.sizes
  S <- dim(all.web.info$predation.matrix)[1]
  ##target.C <- sum(web.matrix[[i]])/length(web.matrix[[i]][,1])^2
  
  for(j in 1:NITS){
    
    ##print(c(i,j))
    
    n = rnorm(1, m.n, sd.n)
    n = 10^n
    ni = rnorm(1, m.ni, sd.ni)
    
    h = rnorm(1, m.h, sd.h)
    h = 10^h
    hi =  rnorm(1, m.hi, sd.hi)
    hj = rnorm(1, m.hj, sd.hj)        
    
    ai = rnorm(1, m.ai, sd.ai)
    aj = rnorm(1, m.aj, sd.aj)
    
    if(!optimise.C){
      a = rnorm(1, m.a, sd.a)
      a = 10^a
    }
    
    if(optimise.C){
      parms <- list(e, n, ni, ai, aj, h, hi, hj, M, S, target.C, real.web)
      a <- get.power.a(parms)
    }
    
    
    EHL <- Allometric.EHL(M,
                          e,
                          h, hi, hj,
                          a, ai, aj,
                          n, ni)
    
    webs[[j]] <- Get.web(EHL)
    EHLs[[j]] <- EHL
    best[[j]] <- Compare.links(all.web.info$predation.matrix, webs[[j]])
  }
  
  list(power=best, EHL=EHLs, web=webs)
}


fit.power.C = function(a, x){
  
  e <- x[[1]]
  n <- x[[2]]
  ni <- x[[3]]
  ai = x[[4]]
  aj = x[[5]]
  h = x[[6]]
  hi = x[[7]]
  hj = x[[8]]
  M = x[[9]]
  S=x[[10]]
  target.C <- x[[11]]
  real.web <- x[[12]]
  
  a <- 10^a
  
  EHL <- Allometric.EHL(M,
                        e,
                        h, hi, hj,
                        a, ai, aj,
                        n, ni)
  web <- Get.web(EHL)
  C=sum(web)/S^2                                             
  ans.wer = abs(target.C - C)
  #print(C)
  ans.wer
}

get.power.C <- function(a.vector, x){
  
  e <- x[[1]]
  n <- x[[2]]
  ni <- x[[3]]
  ai = x[[4]]
  aj = x[[5]]
  h = x[[6]]
  hi = x[[7]]
  hj = x[[8]]
  M = x[[9]]
  S=x[[10]]
  target.C <- x[[11]]
  real.web <- x[[12]]
  
  C = numeric(length=length(a.vector))
  for(i in 1:length(a.vector)){
    a = 10^a.vector[i]
    
    EHL <- Allometric.EHL(M,
                          e,
                          h, hi, hj,
                          a, ai, aj,
                          n, ni)
    
    if(!any(apply(EHL[[1]]/EHL[[2]], 2, var)==0)){
      
      ##      if(!sum(duplicated(apply(EHL[[1]]/EHL[[2]], 2, var)))==(length(EHL[[1]])-1)){
      web <- Get.web(EHL)
      C[i] = sum(web)/S^2
    }
    else{
      web = NA
      C[i] = -1
    }
  }
  C
}    

get.power.a <- function(x){    
  
  e <- x[[1]]
  n <- x[[2]]
  ni <- x[[3]]
  ai = x[[4]]
  aj = x[[5]]
  h = x[[6]]
  hi = x[[7]]
  hj = x[[8]]
  M = x[[9]]
  S=x[[10]]
  target.C <- x[[11]]
  real.web <- x[[12]]
  
  last.a <- 0
  interval <- 0.5
  flag=F
  cant.do=F
  
  ##first check if web has minimum connectance
  if(target.C<(1/S))
    a = 10^last.a
  if(!target.C<(1/S)){
    
    while(flag==F & cant.do==F){
      last.C <- get.power.C(last.a, x)
      ##print(c(last.a, last.C))
      if(last.C==-1){
        flag=T
        cant.do=T
      }
      if(last.C<target.C)
        next.a <- last.a-interval
      if(last.C>target.C)
        next.a <- last.a+interval
      if(last.C==target.C){
        next.a <- last.a-1
        last.a <- last.a+1
        flag=T
      }   
      next.C <- get.power.C(next.a, x)
      if(last.a<next.a & last.C>target.C & next.C<target.C)
        flag=T
      if(last.a>next.a & last.C<target.C & next.C>target.C)
        flag=T
      if(flag==F)
        last.a = next.a
    }
    if(cant.do==F)
      a = 10^optimise(fit.power.C, range(c(last.a, next.a)), x=x)[[1]]
    if(cant.do==T)
      a = NA
  }
  a
}


power.power <- function(opt, x, optimizer, opt.these=c("ai", "aj", "hi", "hj")){
  
  poss.opt <- c("ai", "aj", "hi", "hj")
  
  opt.pars <- c(0,0,0,0)
  
  for(i in 1:length(opt.these))
    opt.pars[match(opt.these[i], poss.opt)] <- opt[i]
  
  ai = opt.pars[1]
  aj = opt.pars[2]
  hi = opt.pars[3]
  hj = opt.pars[4]
  
  
  
  e=x[[1]]
  n=x[[2]]
  ni=x[[3]]
  h=x[[4]]
  M=x[[5]]
  S=x[[6]]
  target.C <- x[[7]]
  real.web <- x[[8]]
  
  parms <- list(e, n, ni, ai, aj, h, hi, hj, M, S, target.C, real.web)
  
  a <- get.power.a(parms)
  
  if(!is.na(a)){
    ##print(c(a, ai, aj, hi, hj))
    
    EHL <- Allometric.EHL(M,
                          e,
                          h, hi, hj,
                          a, ai, aj,
                          n, ni)
    web <- Get.web(EHL)
    C=sum(web)/S^2                                             
    if(abs(target.C-C)>0.05)
      power.yy = 0
    else{
      
      power.yy <- Compare.links(real.web, web)
    }
    if(optimizer)
      result = power.yy
    if(!optimizer)
      result = c(power.yy, a)
  }
  if(is.na(a)){
    if(optimizer)
      result = 0
    if(!optimizer)
      result = c(0, a)
  }
  ##print(result)
  result
}





CE.power <- function(all.web.info){
  
  real.web <- all.web.info$predation.matrix
  
  M <- all.web.info$species.sizes
  S <- dim(real.web)[1]
  target.C <- sum(real.web)/S^2    
  
  e <- 1
  
  n <- 1
  ni <- -3/4
  
  ais <- seq(-1, 1, 0.25)
  ajs <- seq(-1, 1, 0.25)
  
  h <- 1
  his <- seq(1, 2, 0.1)
  hjs <- seq(-2, 0, 0.1)
  
  dont.use <- his!=1
  his <- his[dont.use]    
  
  best = 0
  for(hi in his){
    for(hj in hjs){
      for(ai in ais){
        for(aj in ajs){
          
          ##print(c(hi, hj, ai, aj))
          
          x <- list(e, n, ni, h, M, S, target.C, real.web)
          opt = c(ai, aj, hi, hj)
          
          temp <- power.power(opt=opt, x=x, optimizer=F)
          power.yy <- temp[1]
          a = temp[2]
          if(!is.na(a)){
            EHL <- Allometric.EHL(M,
                                  e,
                                  h, hi, hj,
                                  a, ai, aj,
                                  n, ni)
            web <- Get.web(EHL)
            pred.C = sum(web)/S^2                                             
            
            if(power.yy>best){
              best = power.yy
              power.results = c(target.C, pred.C, hi, hj, ai, aj, power.yy, a)
              best.web = web
              best.EHL = EHL
            }
          }
        }
      }
    }
  }
  
  list(power=best, pars=c(n=n, ni=ni, h=h, hi=hi, hj=hj, a=a, ai=ai, aj=aj), EHL=EHL, web=best.web)
}




NM.power <- function(all.web.info){
  
  if(.Platform$OS.type=="windows")
    power.initial.pars <- as.matrix(read.csv("C:\\Documents and Settings\\bo1op\\My Documents\\work\\research\\3.in.progress\\allometric.web\\data\\NM.optim\\power.initial.pars.txt", header=F))
  if(.Platform$OS.type=="unix")
    power.initial.pars <- as.matrix(read.csv("~/allometric.web/data/NM.optim/power.initial.pars.txt", header=T))
  
  real.web <- all.web.info$predation.matrix
  M <- all.web.info$species.sizes
  S <- dim(real.web)[1]
  target.C <- sum(real.web)/S^2    
  
  e <- 1
  n <- 1
  ni <- -3/4
  h <- 1
  perm.parms <- list(e=e, n=n, ni=ni, h=h, M=M, S=S, target.C=target.C, real.web)
  
  best = 0
  
  for(ip in 1:length(power.initial.pars[,1])){
    
    power.power(opt=power.initial.pars[ip,], x=perm.parms, optimizer=T)
    
    o.p <- optim(power.initial.pars[ip,], power.power, control=list(fnscale=-1, trace=0), x=perm.parms, optimizer=T)
    print(o.p)
    parms <- list(e, n, ni, ai=o.p$par[1], aj=o.p$par[2], h, hi=o.p$par[3], hj=o.p$par[4], M, S=S, target.C=target.C, real.web)
    a = get.power.a(parms)
    if(o.p$value>best){
      best=o.p$value
      optim.power.pars <- c(e=e, n=n, ni=ni, a=a, ai=o.p$par[1], aj=o.p$par[2], h=h, hi=o.p$par[3], hj=o.p$par[4])           
    }
  }
  
  EHL <- Allometric.EHL(M,
                        e,
                        h, hi=optim.power.pars[8], hj=optim.power.pars[9],
                        a=optim.power.pars[4], ai=optim.power.pars[5], aj=optim.power.pars[6], 
                        n, ni)
  web <- Get.web(EHL)
  print(sum(web))
  
  list(power=best, pars=optim.power.pars, EHL=EHL, web=web)
  write.table(optim.power.pars, file=paste("power.optim.pars.txt", sep="."))
  
}

NM.power.comp <- function(all.web.info){
  
  if(.Platform$OS.type=="windows")
    power.initial.pars <- as.matrix(read.csv("C:\\Documents and Settings\\bo1op\\My Documents\\work\\research\\3.in.progress\\allometric.web\\data\\NM.optim\\power.initial.pars.txt", header=F))
  if(.Platform$OS.type=="unix")
    power.initial.pars <- as.matrix(read.csv("~/allometric.web/data/NM.optim/power.initial.pars.txt", header=T))
  
  real.web <- all.web.info$predation.matrix
  M <- all.web.info$species.sizes
  S <- dim(real.web)[1]
  target.C <- sum(real.web)/S^2    
  
  e <- 1
  n <- 1
  ni <- -3/4
  h <- 1
  perm.parms <- list(e=e, n=n, ni=ni, h=h, M=M, S=S, target.C=target.C, real.web)
  
  best = 0
  
  opt.pars <- c("ai", "aj", "hi", "hj")
  
  par.holds <- cbind(rep(c(0,1), each=8, length.out=16),
                     rep(c(0,1), each=4, length.out=16),
                     rep(c(0,1), each=2, length.out=16),
                     rep(c(0,1), each=1, length.out=16))
  
  best.by.ph <- list()
  
  for(ph in 2:length(par.holds[,1])){
    print(ph)
    best = -0.1
    
    opt.these <- opt.pars[par.holds[ph,]==1]
    
    for(ip in 1:length(power.initial.pars[,1])){
      
      ##power.power(opt=power.initial.pars[ip,], x=perm.parms, optimizer=T)
      
      print(c(ph, ip))
      
      if(length(opt.these)>1)
        o.p <- optim(power.initial.pars[ip,par.holds[ph,]==1],
                     power.power,
                     control=list(fnscale=-1, trace=1),
                     x=perm.parms, optimizer=T,
                     opt.these=opt.these)
      
      if(length(opt.these)==1){
        o.p <- list(par=NA, value=NA)
        opt.1d <- optimise(power.power,
                           lower=-2, upper=2,
                           maximum=T,
                           x=perm.parms, optimizer=T,
                           opt.these=opt.these)
        o.p$par <- opt.1d$maximum
        o.p$value <- opt.1d$objective
      }            
      
      
      poss.opt <- opt.pars
      
      opt1.pars <- c(0,0,0,0)
      
      for(i in 1:length(opt.these))
        opt1.pars[match(opt.these[i], poss.opt)] <- o.p$par[i]
      
      ai = opt1.pars[1]
      aj = opt1.pars[2]
      hi = opt1.pars[3]
      hj = opt1.pars[4]
      
      parms <- list(e, n, ni,
                    ai=ai, aj=aj,
                    h, hi=hi, hj=hj,
                    M, S=S, target.C=target.C, real.web)
      
      a = get.power.a(parms)
      
      if(o.p$value>best){
        best=o.p$value
        optim.power.pars <- c(e=e, n=n, ni=ni,
                              a=a, ai=ai, aj=aj,
                              h=h, hi=hi, hj=hj)           
      }
    }        
    EHL <- Allometric.EHL(M,
                          e,
                          h, hi=optim.power.pars[8], hj=optim.power.pars[9],
                          a=optim.power.pars[4], ai=optim.power.pars[5], aj=optim.power.pars[6], 
                          n, ni)
    web <- Get.web(EHL)
    ##print(sum(web))
    
    ##list(power=best, pars=optim.power.pars, EHL=EHL, web=web)
    ##write.table(optim.power.pars, file=paste("power.optim.pars.txt", sep="."))
    
    best.by.ph[[ph]] <- list(best, par.holds[ph,], web)
    
  }
  
  best.by.ph
}


fit.ratio.C = function(a, x){
  
  e <- x[[1]]
  n <- x[[2]]
  ni <- x[[3]]
  ai = x[[4]]
  aj = x[[5]]
  r.a = x[[6]]
  r.b = x[[7]]
  M = x[[8]]
  S=x[[9]]
  target.C <- x[[10]]
  real.web <- x[[11]]
  
  
  a <- 10^a
  
  EHL <- Ratio.allometric.EHL(M=M,
                              e=e,
                              a=a, ai=ai, aj=aj,
                              n=n, ni=ni,
                              r.a=r.a, r.b=r.b)
  web <- Get.web(EHL)
  C=sum(web)/S^2                                             
  ans.wer = abs(target.C - C)
  #print(C)
  ans.wer
}

get.ratio.C <- function(a.vector, x){
  
  e <- x[[1]]
  n <- x[[2]]
  ni <- x[[3]]
  ai = x[[4]]
  aj = x[[5]]
  r.a = x[[6]]
  r.b = x[[7]]
  M = x[[8]]
  S=x[[9]]
  target.C <- x[[10]]
  real.web <- x[[11]]
  
  C = numeric(length=length(a.vector))
  for(i in 1:length(a.vector)){
    a = 10^a.vector[i]
    
    EHL <- Ratio.allometric.EHL(M=M,
                                e=e,
                                a=a, ai=ai, aj=aj,
                                n=n, ni=ni,
                                r.a=r.a, r.b=r.b)
    
    web <- Get.web(EHL)
    C[i] = sum(web)/S^2                                             
  }
  C
}    


get.ratio.a <- function(x){
  
  e <- x[[1]]
  n <- x[[2]]
  ni <- x[[3]]
  ai = x[[4]]
  aj = x[[5]]
  r.a = x[[6]]
  r.b = x[[7]]
  M = x[[8]]
  S=x[[9]]
  target.C <- x[[10]]
  real.web <- x[[11]]
  
  last.a <- 0
  interval <- 0.5
  flag=F
  
  ## first check if web has minimum connectance
  if(target.C<(1/S))
    a = 10^last.a
  if(!target.C<(1/S)){
    
    while(flag==F){
      last.C <- get.ratio.C(last.a, x)
      if(last.C<target.C)
        next.a <- last.a-interval
      if(last.C>target.C)
        next.a <- last.a+interval
      if(last.C==target.C){
        next.a <- last.a-1
        last.a <- last.a+1
        flag=T
      }   
      next.C <- get.ratio.C(next.a, x)
      if(last.a<next.a & last.C>target.C & next.C<target.C)
        flag=T
      if(last.a>next.a & last.C<target.C & next.C>target.C)
        flag=T
      if(flag==F)
        last.a = next.a
    }    
    a = 10^optimise(fit.ratio.C, range(c(last.a, next.a)), x=x)[[1]]
  }
  a
}
