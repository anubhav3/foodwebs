#rm(list=ls())
## This function takes the model parameters and 
## body masses, and returns the
## vector of energies E
## matrix of handling times H
## and matrix of encounter rates L
Ratio.allometric.EHL <- function(M,
                                 e,
                                 r.a, r.b,
                                 a, ai, aj,
                                 n, ni=-3/4){
  
  ## The handling time function
  get.h <- function(Mi, Mj, r.a, r.b)
    ifelse((r.b-Mi/Mj)>0,
           r.a/(r.b-Mi/Mj),
           Inf)
  
  
  ## in matrix H resources are rows and consumers are columns
  if(!r.b==0)
    H <- outer(M, M, get.h, r.a, r.b)
  if(r.b==0)
    H = matrix(r.a, length(M),length(M))
  
  
  ## ENCOUNTER RATES: consumer - resource mass specific encounter rates
  get.a <- function(Mi, Mj,
                    a, ai, aj)
    a * Mi^ai * Mj^aj
  A <- outer(M, M, get.a,
             a=a, ai=ai, aj=aj)
  L <- A* n*M^ni
  
  ## Check if body sizes are sorted, as they need to be
  if(sum(order(M)==1:length(M))!=length(M))
    stop("Body sizes not sorted")
  
  ## energy values
  E <- e*M
  
  ## object to return
  list(E=E, H=H, L=L)
  
}


## Function taking the energy, handling time, encounter rate object
## and returns the binary food web
## and total energy intake rate of speies as consumers if energy.intake == TRUE
Get.web <- function(EHL, energy.intake=F){
  
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
      
      p <- P[,j]
      
      if(sum(p>0)==1)
        web[which(p>0),j] <- 1
      
      if(sum(p>0)>1){
        
        ## ordering of p required
        
        order.by.p <- order(p, decreasing=T)
        p <- p[order.by.p]
        Lj <- EHL[[3]][,j][order.by.p]
        hj <- EHL[[2]][,j][order.by.p]
        Ej <- EHL[[1]][order.by.p]
        
        cumulative.profit <- cumsum(Ej * Lj) / (1 + cumsum( Lj * hj))
        
        dj <- max(which(cumulative.profit==max(cumulative.profit)))
        
        web[,j] <- c(rep(1, dj), rep(0, S-dj))[order(order.by.p)]
        
        overall.energy[j] <- cumulative.profit[dj]
        
        energies <- c(Ej * Lj)[1:sum(web[,j])] / (1 + cumsum( Lj * hj)[sum(web[,j])]) 
        all.energies <- c(energies, rep(0, S-length(energies)))
        
        per.species.energy[,j] <- all.energies[order(order.by.p)]
        
      }
    }
  }
  
  ## This if encounter rates are not predator specific
  if(is.vector(EHL[[3]])){
    for(j in 1:S){
      
      if(sum(p>0)==1)
        web[which(p>0),j] <- 1
      
      if(sum(p>0)>1){
        
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
        
        overall.energy[j] <- cumulative.profit[dj]    
        
        energies <- c(Ej * Lj)[1:sum(web[,j])] / (1 + cumsum( Lj * hj)[sum(web[,j])]) 
        all.energies <- c(energies, rep(0, S-length(energies)))
        
        per.species.energy[,j] <- all.energies[order(order.by.p)]
        
      }
    }
  }
  
  if(energy.intake)
    result <- list(web=web, overall.flux=overall.energy, per.species.flux=per.species.energy)
  else
    result <- web
  result
  
}


## the function that takes the parameters to optimise in opt and
## returns the number of links
ratio.power <- function(opt, x){
  
  #print(opt)
  
  ## recall that parameters a and r.b are estimated on the exponential scale
  ## so the real value is 10^the estimate value.
  
  a = 10^as.numeric(opt["a"])
  ai = as.numeric(opt["ai"])
  aj = as.numeric(opt["aj"])
  r.b = 10^as.numeric(opt["r.b"])
  
  e=x[["e"]]
  n=x[["n"]]
  ni=x[["ni"]]
  r.a=x[["r.a"]]
  M=x[["M"]]
  S=x[["S"]]
  target.C <- x[["target.C"]]
  real.web <- x[["real.web"]]
  
  #parms <- list(e, n, ni, ai, aj, r.a, r.b, M, S, target.C, real.web)
  
  EHL <- Ratio.allometric.EHL(M=M,
                              e=e,
                              a=a, ai=ai, aj=aj,
                              n=n, ni=ni,
                              r.a=r.a, r.b=r.b)
  
  ##print(sum(EHL[[2]]!=Inf))
  
  # if(sum(EHL[[2]]!=Inf)>=(target.C*S^2)){
  #   a = get.ratio.a(parms)
  #   EHL <- Ratio.allometric.EHL(M=M,
  #                               e=e,
  #                               a=a, ai=ai, aj=aj,
  #                               n=n, ni=ni,
  #                               r.a=r.a, r.b=r.b)
  web <- Get.web(EHL)
  #  C=sum(web)/S^2                                             
  #  if(abs(target.C-C)>0.05)
  #    ratio.yy = -1
  #  else{            
  #ratio.yy <- Compare.links(real.web, web)
  #}
  #}
  #if(sum(EHL[[2]]!=Inf)<(target.C*S^2))
  #  ratio.yy = -1    
  #if(optimizer)
    result = sum(web)  #It will return the number of links
  #if(!optimizer)
  #  result = c(ratio.yy, a)
  ##print(paste("XXX", result))
  result
}




## The function to call in order to do the optimisation
abc.ratio <- function(all.web.info, starting_parameter_values, tol=tol){
  
  ## the function runs the optimisation separately for each row of
  ## parameter starting values
  ratio.initial.pars <- starting_parameter_values
  
  real.web <- all.web.info$predation.matrix
  M <- all.web.info$species.sizes
  S <- dim(real.web)[1]
  target.C <- sum(real.web)/S^2    
  e <- 1
  n <- 1
  ni <- -3/4
  r.a = 1
  parms <- list(e=e, n=n, ni=ni, r.a=r.a, M=M, S=S, target.C=target.C, real.web=real.web)
  ##ratio.power(ratio.initial.pars[ip,], parms)
  
  
  n.sample = length(ratio.initial.pars[,1])
  dist = rep(NA,n.sample)
  for(i in 1:n.sample){
    dist[i] = abs(ratio.power(opt=ratio.initial.pars[i,], x=parms)/S^2-sum(real.web)/S^2)
  } 
  
  post_dists <- ratio.initial.pars[dist<tol,]
  
  pred.means = colMeans(post_dists)
  pred.sd = apply(post_dists,2,sd)
  
  
  best.EHL <- Ratio.allometric.EHL(M=M,
                                   e=e,
                                   a=10^pred.means[1], ai=pred.means[2], aj=pred.means[3],
                                   n=n, ni=ni,
                                   r.a=r.a, r.b=10^pred.means[4])
  best.web <- Get.web(best.EHL) 
  optim.power.pars <- c(e=e, n=n, ni=ni, a=10^pred.means[1], ai=pred.means[2],
                        aj=pred.means[3], r.a=r.a, r.b=10^pred.means[4])
  acc = ratio.power(opt = pred.means , x=parms)
  
  list(acc = acc ,post = post_dists, pars=optim.power.pars, EHL=best.EHL, web=best.web)
}



Compare.links <- function(real.web, model.web){
  
  ## Here is the original function used in the PNAS 2008 paper
  ## result <- sum(real.web==1 & model.web==1)  / sum(model.web)
  
  ## New function, accuracy
  result <- sum(real.web == model.web) / dim(real.web)[[1]]^2
  
  result
}
