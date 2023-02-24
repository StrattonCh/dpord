##############
## Fit CORAL models using JAGS
## Usually d = 0-2 latent variables. The number of clusters, K, can take any non-negative integer
## Binomial, Negative Binomial, and Poisson distributions allowed
## NOTE: Appropriate adjustments needs to be made to the JAGS scripts in coral-auxilaryfunctions.R to handle CORAL models with d == 1 and/or K == 1
##############
library(R2jags); 
library(mvtnorm); 
library(label.switching); 
library(gtools); 
library(coda)
source("coral-auxilaryfunctions.R")

fit.coral <- function(y, family, num.lv = 2, num.clus = 2, site.eff = TRUE, n.burnin = 5000, n.iter = 55000, n.thin = 10, seed = 123, labswitch.method = "ecr", init.mod = NULL, save.model = FALSE) {  
  if(num.lv >= ncol(y)) stop("# of LVs >= # of y...doesn't make sense.")
  if(num.clus == 1) print("CORAL used for pure ordination.")
  if(num.lv == 0) print("No latent variables included; num.clus argument ignored.")
  if(!is.null(init.mod)) { print("Initial fit provided.") }
  
  site.eff <- as.numeric(site.eff)
  N <- nrow(y); S <- ncol(y)
  
  ## If d == 0, then make things easy simple
  if(num.lv == 0) {
    jags.data <- list(y = y, N = N, S = S, num.lv = num.lv, site.eff = site.eff)
    jags.params <- c("all.params")
    if(site.eff) jags.params <- c(jags.params, "site.params")
    set.seed(seed)
    if(family == "binomial")
      jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.binom.lv0clus0, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) 
    if(family == "poisson")
      jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.pois.lv0clus0, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) 
    if(family == "negative.binomial")
      jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.nb.lv0clus0, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) 
    
    ## Format into array
    fit.mcmcBase <- jagsfit$BUGSoutput
    jagsfit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = n.thin) 
    combined.jagsfit.mcmc <- jagsfit.mcmc
    rm(jagsfit.mcmc)
    
    fit <- list(
      params.med = matrix(apply(combined.jagsfit.mcmc[,grep("all.params",colnames(combined.jagsfit.mcmc))],2,median),nrow=S),
      params.mean = matrix(apply(combined.jagsfit.mcmc[,grep("all.params",colnames(combined.jagsfit.mcmc))],2,mean),nrow=S),
      site.params.med = site.params.mean <- NULL)
    
    if(site.eff) {
      fit$site.params.med <- apply(combined.jagsfit.mcmc[,grep("site.params",colnames(combined.jagsfit.mcmc))],2,median)
      fit$site.params.mean <- apply(combined.jagsfit.mcmc[,grep("site.params",colnames(combined.jagsfit.mcmc))],2,mean) }
    
    if(save.model) fit$jags.model <- jagsfit
    
    fit$num.lv <- num.lv; 
    fit$num.clus <- num.clus
    fit$y <- y; 
    fit$family <- family
    fit$site.eff <- site.eff
    fit$n.burnin <- n.burnin;
    fit$n.iter <- n.iter
    fit$n.thin <- n.thin
    
    return(fit)
  }
  
  ## Fit an inital model to obtain sites to pivot CORAL around -- necessary for JAGS/BUGS
  if(is.null(init.mod)) {
    set.seed(seed); try.seed.counter <- 0; try.counter <- 0
    init.mod <- try(fit.initmod(y = y, family = family, num.lv = num.lv, site.eff = site.eff, n.burnin = 5000, n.iter = 20000, n.thin = 5),silent=T)
    
    while(inherits(init.mod,"try-error") & try.counter < 5) {
      try.seed.counter <- try.seed.counter + n.burnin
      set.seed(try.seed.counter)
      init.mod <- try(fit.initmod(y = y, family = family, num.lv = num.lv, site.eff = site.eff, n.burnin = 5000, n.iter = 15000, n.thin = 5),silent=T)
      try.counter <- try.counter + 1
    }
    
    if(try.counter >= 5)  { cat("Failure to initialize model...exiting. Sorry!"); break; } 
  }
  
  ## If K == 1, then stop and done!
  if(num.clus == 1) {
    fit <- init.mod 
    fit$mu.clus.med <- fit$mu.clus.mean <- matrix(0,1,num.lv)
    fit$pis.med <- fit$pis.mean <- 1
    
    fit$num.lv <- num.lv; 
    fit$num.clus <- num.clus
    fit$y <- y; 
    fit$family <- family
    fit$site.eff <- site.eff
    fit$n.burnin <- n.burnin;
    fit$n.iter <- n.iter
    fit$n.thin <- n.thin
    
    return(fit)
  }
  
  cluster.lvs <- kmeans(init.mod$lv.med, centers = num.clus, nstart = 100)
  find.closest.points <- matrix(NA,num.clus,2)
  for(k in 1:num.clus) { 
    calc.dist <- as.matrix(dist(rbind(cluster.lvs$centers[k,],init.mod$lv.med)))[-1,1] ## Calculate euclidian distances from cluster center to all points 
    find.closest.points[k,] <- order(calc.dist)[1:2] } ## Find the two sites with minimum distance to cluster mean, and use them as pivot sites
  make.clus.labels <- rep(NA,N); 
  make.clus.labels[as.vector(find.closest.points)] <- rep(1:num.clus,2) 
  
  jags.data <- list(y = y, N = N, S = S, num.lv = num.lv, num.clus = num.clus, site.eff = site.eff, alpha = rep(1,num.clus), clus.label = make.clus.labels)
  jags.params <- c("all.params","lvs","clus.label","mixprop","mu.clus")
  if(num.clus == 1) { 
    jags.data <- list(y = y, N = N, S = S, num.lv = num.lv, site.eff = site.eff)
    jags.params <- c("all.params","lvs","mu.clus") }
  if(site.eff) { jags.params <- c(jags.params, "site.params") }
  set.seed(seed)
  
  ## Fit the model
  if(family == "binomial") {
    if(num.lv == 1 & num.clus > 1) jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.binom.lv1, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) 
    if(num.lv > 1 & num.clus > 1) jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.binom.lv, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) }
  if(family == "poisson") {
    if(num.lv == 1 & num.clus > 1) jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.pois.lv1, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) 
    if(num.lv > 1 & num.clus > 1) jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.pois.lv, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) }
  if(family == "negative.binomial") {
    if(num.lv == 1 & num.clus > 1) jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.nb.lv1, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) 
    if(num.lv > 1 & num.clus > 1) jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = coralmod.nb.lv, n.iter = n.iter, n.burnin = n.burnin, n.thin = n.thin, n.chains = 1, DIC = T) }
  
  ## Format into array
  fit.mcmcBase <- jagsfit$BUGSoutput
  jagsfit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = n.thin) 
  combined.jagsfit.mcmc <- jagsfit.mcmc
  rm(jagsfit.mcmc)
  
  ## Apply Equivalence Classes algorithm
  if(labswitch.method == "ecr" & num.clus > 1) {
    cat("ECR Label switching algorithm applied.\n")
    lv.med <- matrix(apply(combined.jagsfit.mcmc[,grep("lvs",colnames(combined.jagsfit.mcmc))],2,mean),nrow=N)
    get.kmeans <- kmeans(lv.med, centers = num.clus, nstart = 100)
    
    params.array <- array(NA,dim=c(nrow(combined.jagsfit.mcmc), num.clus, num.lv+1)) ## + 1 for mix prop
    params.array[,,1] <- combined.jagsfit.mcmc[,grep("mixprop",colnames(combined.jagsfit.mcmc))]
    params.array[,,2:(num.lv+1)] <- combined.jagsfit.mcmc[,grep("mu.clus",colnames(combined.jagsfit.mcmc))]
    
    run.ecr <- ecr(zpivot = get.kmeans$clus, z = combined.jagsfit.mcmc[,grep("clus.label",colnames(combined.jagsfit.mcmc))], K = num.clus)
    sub.jagsfit.mcmc <- permute.mcmc(params.array, run.ecr$permutations)$output
    all.z <- permuted.z <- combined.jagsfit.mcmc[,grep("clus.label",colnames(combined.jagsfit.mcmc))]
    for(k in 1:nrow(combined.jagsfit.mcmc)) { for(l in 1:num.clus) { permuted.z[k,all.z[k,]==l] = run.ecr$permutations[k,l] } }
    rm(lv.med, get.kmeans)
  }
  
  ## Apply Stephens' reordering algorithm
  if(labswitch.method == "stephens" & num.clus > 1) {
    cat("Stephens Label switching algorithm applied\n")
    all.taus <- array(NA, dim = c(nrow(combined.jagsfit.mcmc), N, num.clus))
    for(t in 1:nrow(combined.jagsfit.mcmc)) { 
      all.taus[t,,] <- calc.taus(data = y, mixprops = combined.jagsfit.mcmc[t,grep("mixprop",colnames(combined.jagsfit.mcmc))], 
                                 mu.clus = matrix(combined.jagsfit.mcmc[t,grep("mu.clus",colnames(combined.jagsfit.mcmc))],num.clus), 
                                 lvs = matrix(combined.jagsfit.mcmc[t,grep("lvs",colnames(combined.jagsfit.mcmc))],N))$taus }
    
    params.array <- array(NA,dim=c(nrow(combined.jagsfit.mcmc), num.clus, num.lv+1)) ## + 1 for mix prop
    params.array[,,1] <- combined.jagsfit.mcmc[,grep("mixprop",colnames(combined.jagsfit.mcmc))]
    params.array[,,2:(num.lv+1)] <- combined.jagsfit.mcmc[,grep("mu.clus",colnames(combined.jagsfit.mcmc))]
    
    run.stephens <- stephens(p = all.taus)
    sub.jagsfit.mcmc <- permute.mcmc(params.array, run.stephens$permutations)$output
    all.z <- permuted.z <- combined.jagsfit.mcmc[,grep("clus.label",colnames(combined.jagsfit.mcmc))]
    for(k in 1:nrow(combined.jagsfit.mcmc)) { for(l in 1:num.clus) { permuted.z[k,all.z[k,]==l] = run.stephens$permutations[k,l] } }
    rm(all.taus)
  }
  
  ## Apply no label switching
  if((labswitch.method == "none" & num.clus > 1) | num.clus == 1) {
    cat("No Label switching algorithm applied\n")
    sub.jagsfit.mcmc <- array(NA,dim=c(nrow(combined.jagsfit.mcmc), num.clus, num.lv+1)) ## + 1 for mix prop
    if(num.clus > 1) sub.jagsfit.mcmc[,,1] <- combined.jagsfit.mcmc[,grep("mixprop",colnames(combined.jagsfit.mcmc))]
    if(num.clus == 1) sub.jagsfit.mcmc[,,1] <- 1
    sub.jagsfit.mcmc[,,2:(num.lv+1)] <- combined.jagsfit.mcmc[,grep("mu.clus",colnames(combined.jagsfit.mcmc))]	
    if(num.clus > 1) permuted.z <- combined.jagsfit.mcmc[,grep("clus.label",colnames(combined.jagsfit.mcmc))] 
  }
  
  mu.clus.med <- mu.clus.mean <- matrix(NA,num.clus,num.lv) 
  for(i in 1:num.clus) { for(j in 1:num.lv) { 
    mu.clus.med[i,j] <- median(sub.jagsfit.mcmc[,i,j+1]); mu.clus.mean[i,j] <- mean(sub.jagsfit.mcmc[,i,j+1]) } }
  
  fit <- list(
    mu.clus.med = mu.clus.med,
    params.med = matrix(apply(combined.jagsfit.mcmc[,grep("all.params",colnames(combined.jagsfit.mcmc))],2,median),nrow=S),
    lv.med = matrix(apply(combined.jagsfit.mcmc[,grep("lvs",colnames(combined.jagsfit.mcmc))],2,median),nrow=N), 
    
    mu.clus.mean = mu.clus.mean,
    params.mean = matrix(apply(combined.jagsfit.mcmc[,grep("all.params",colnames(combined.jagsfit.mcmc))],2,mean),nrow=S),
    lv.mean = matrix(apply(combined.jagsfit.mcmc[,grep("lvs",colnames(combined.jagsfit.mcmc))],2,mean),nrow=N)
  )
  
  fit$pis.med = apply(sub.jagsfit.mcmc[,,1],2,median); fit$pis.mean = apply(sub.jagsfit.mcmc[,,1],2,mean)
  fit$clus.modal = apply(permuted.z,2,get.mode); fit$clus.table = apply(permuted.z,2,table)
  
  if(site.eff) {
    fit$site.params.med = apply(combined.jagsfit.mcmc[,grep("site.params",colnames(combined.jagsfit.mcmc))],2,median)
    fit$site.params.mean = apply(combined.jagsfit.mcmc[,grep("site.params",colnames(combined.jagsfit.mcmc))],2,mean) }
  if(!site.eff) { fit$site.params.med <- fit$site.params.mean <- NULL }
  
  if(save.model) fit$jags.model <- jagsfit
  
  fit$init.fit <- init.mod
  fit$find.closest.points <- find.closest.points
  fit$num.lv <- num.lv; 
  fit$num.clus <- num.clus
  fit$y <- y; 
  fit$family <- family
  fit$site.eff <- site.eff
  fit$labswitch.method <- labswitch.method
  fit$n.burnin <- n.burnin;
  fit$n.iter <- n.iter
  fit$n.thin <- n.thin
  
  return(fit) 
}


get.ics <- function(coralfit) {
  y <- coralfit$y
  family <- coralfit$family
  n <- nrow(y); s <- ncol(y)
  
  if(coralfit$num.lv == 0) {
    print("No latent variables detected in coralfit")
    do.ics <- get.ics.lv0(coralfit = coralfit)
    return(do.ics)	
  }
  
  if(is.null(coralfit$mu.clus.mean)) { 
    num.lv <- ncol(coralfit$lv.med)
    coralfit$pis.med <- coralfit$pis.mean <- 1
    coralfit$mu.clus.med <- coralfit$mu.clus.mean <- matrix(0,1,num.lv) }
  
  num.lv <- ncol(coralfit$mu.clus.mean); num.clus <- nrow(coralfit$mu.clus.mean)
  mlogl.med <- get.marg.predlogl(y = y, params = coralfit$params.med, site.params = coralfit$site.params.med, pis = coralfit$pis.med, 
                                 mu.clus = coralfit$mu.clus.med, family)
  clogl.med <- get.cond.predlogl(y = y, lvs = coralfit$lv.med, params = coralfit$params.med, site.params = coralfit$site.params.med, 
                                 family = family, phi = T)$predlogl.mat
  
  num.params <- s*num.lv - 0.5*num.lv*(num.lv-1) + s + s*as.numeric(family == "negative.binomial") + n*(!is.null(coralfit$site.params.mean)) + (num.clus-1)*num.lv + (num.clus-1)
  
  aic.post.med <- -2*sum(mlogl.med[is.finite(mlogl.med)]) + 2*num.params
  bic.post.med <- -2*sum(mlogl.med[is.finite(mlogl.med)]) + log(n)*num.params
  aic.post.condmed <- -2*sum(clogl.med[is.finite(clogl.med)]) + 2*num.params
  bic.post.condmed <- -2*sum(clogl.med[is.finite(clogl.med)]) + log(n)*num.params
  
  if(!is.null(coralfit$jags.model)) {
    get.dic <- coralfit$jags.model$BUGSoutput$DIC
    
    fit.mcmcBase <- coralfit$jags.model$BUGSoutput
    combined.jagsfit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = coralfit$n.thin) 
    n.iterations <- nrow(combined.jagsfit.mcmc)
    all.cond.logl <- numeric(n.iterations)
    
    for(t in 1:n.iterations) {
      cw.params <- matrix(combined.jagsfit.mcmc[t,grep("all.params",colnames(combined.jagsfit.mcmc))],nrow=s)
      cw.lvs <- matrix(combined.jagsfit.mcmc[t,grep("lvs",colnames(combined.jagsfit.mcmc))],nrow=n) 
      if(!is.null(coralfit$site.params.mean)) cw.site.params <- combined.jagsfit.mcmc[t,grep("site.params",colnames(combined.jagsfit.mcmc))]
      if(is.null(coralfit$site.params.mean)) cw.site.params <- NULL
      clogl <- get.cond.predlogl(y = y, lvs = cw.lvs, params = cw.params, site.params = cw.site.params, family = family, phi = T)
      all.cond.logl[t] <- sum(clogl$predlogl.mat) }
    
    eaic <- -2*mean(all.cond.logl[is.finite(all.cond.logl)]) + 2*num.params
    ebic <- -2*mean(all.cond.logl[is.finite(all.cond.logl)]) + log(n)*num.params
  }
  
  ics.out <- c(ifelse(exists("get.dic"), get.dic, NA), aic.post.condmed, bic.post.condmed, aic.post.med, bic.post.med, ifelse(exists("eaic"), eaic, NA), ifelse(exists("ebic"), ebic, NA))
  names(ics.out) <- c("DIC_c","AIC_c,med","BIC_c,med","AIC_med","BIC_med","EAIC_c","EBIC")	
  
  return(list(marg.logl.med = sum(mlogl.med[is.finite(mlogl.med)]), cond.logl.mean = sum(clogl.med[is.finite(clogl.med)]), num.params = num.params, ics = ics.out)) 
}	


## Internal function for calculating logL for d = 0 CORAL models 
## For null models, conditional and marginal log-likelihood are the same
calc.logLik.lv0 <- function(y, family, coefs, site.coefs = NULL) {
  n <- nrow(y); p <- ncol(y); 
  logl <- 0; logl.comp <- numeric(n)
  if(is.null(site.coefs)) site.coefs <- rep(0,n)
  
  for(i in 1:n) {
    eta <- site.coefs[i] + coefs[,1]
    if(family == "poisson") logl.comp[i] <- sum(dpois(as.vector(unlist(y[i,])), lambda = exp(eta), log = T))
    if(family == "binomial") logl.comp[i] <- sum(dbinom(as.vector(unlist(y[i,])), 1, prob = exp(eta)/(1+exp(eta)), log = T))
    if(family == "negative.binomial") logl.comp[i] <- sum(dnbinom(as.vector(unlist(y[i,])), mu = exp(eta), size = 1/coefs[,2], log = T)) }	 	
  
  return(list(logLik = sum(logl.comp), logLik.comp = logl.comp)) }


## Internal function for calculating ics for d = 0 CORAL models 
## Run inside get.ics
get.ics.lv0 <- function(coralfit) {
  y <- coralfit$y
  family <- coralfit$family
  n <- nrow(y); s <- ncol(y)
  
  mlogl.med <- clogl.med <- calc.logLik.lv0(y, family, coefs = coralfit$params.med, site.coefs = coralfit$site.params.med)$logLik.comp
  num.params <- s + s*as.numeric(family == "negative.binomial") + n*(!is.null(coralfit$site.params.mean)) 
  
  aic.post.med <- aic.post.condmed <- -2*sum(mlogl.med[is.finite(mlogl.med)]) + 2*num.params
  bic.post.med <- bic.post.condmed <- -2*sum(mlogl.med[is.finite(mlogl.med)]) + log(n)*num.params
  
  if(!is.null(coralfit$jags.model)) {
    get.dic <- coralfit$jags.model$BUGSoutput$DIC
    
    fit.mcmcBase <- coralfit$jags.model$BUGSoutput
    combined.jagsfit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = coralfit$n.thin) 
    n.iterations <- nrow(combined.jagsfit.mcmc)
    all.cond.logl <- numeric(n.iterations)
    
    for(t in 1:n.iterations) {
      cw.params <- matrix(combined.jagsfit.mcmc[t,grep("all.params",colnames(combined.jagsfit.mcmc))],nrow=s)
      if(!is.null(coralfit$site.params.mean)) cw.site.params <- combined.jagsfit.mcmc[t,grep("site.params",colnames(combined.jagsfit.mcmc))]
      if(is.null(coralfit$site.params.mean)) cw.site.params <- NULL
      clogl <- calc.logLik.lv0(y, family, coefs = cw.params, site.coefs = cw.site.params)$logLik.comp
      all.cond.logl[t] <- sum(clogl) }
    
    eaic <- -2*mean(all.cond.logl[is.finite(all.cond.logl)]) + 2*num.params
    ebic <- -2*mean(all.cond.logl[is.finite(all.cond.logl)]) + log(n)*num.params
  }
  
  ics.out <- c(ifelse(exists("get.dic"), get.dic, NA), aic.post.condmed, bic.post.condmed, aic.post.med, bic.post.med, ifelse(exists("eaic"), eaic, NA), ifelse(exists("ebic"), ebic, NA))
  names(ics.out) <- c("DIC_c","AIC_c,med","BIC_c,med","AIC_med","BIC_med","EAIC_c","EBIC")	
  
  return(list(marg.logl.med = sum(mlogl.med[is.finite(mlogl.med)]), cond.logl.mean = sum(mlogl.med[is.finite(mlogl.med)]), 
              num.params = num.params, ics = ics.out)) 
}	
