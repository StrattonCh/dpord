###############
## K = 1 CORAL models, fitted also to obtain "pivot" sites - LVs drawn from independent standard gaussian
###############
initmod.nb.lv <- function() {
  ## Likelihood
  for(j in 1:s) { for(i in 1:n) { 
    eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
    mu[i,j] <- exp(eta[i,j])
    u[i,j] ~ dgamma(all.params[j,num.lv+2], all.params[j,num.lv+2])
    y[i,j] ~ dpois(mu[i,j]*u[i,j]) ## Parameterizing the NB as a multiplicative random effect models
  } }
  
  for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } }
  
  ## Prior
  for(j in 1:s) { all.params[j,1] ~ dnorm(0,1e-2) } ## species-specific intercept	
  for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { all.params[i,j] <- 0 } } ## Upper diagonal constrained to 0
  for(i in 1:num.lv) { all.params[i,i+1] ~ dunif(0,1e3) } ## Constraints on diagonal 	
  for(i in 2:num.lv) { for(j in 2:i) { all.params[i,j] ~ dnorm(0,1e-2) } } ## Lower diagonal unconstrained
  for(i in (num.lv+1):s) { for(j in 2:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } ## All other elements
  for(j in 1:s) { all.params[j,num.lv+2] ~ dunif(0,1e3) }
  
  for(i in 1:n) { site.params[i] ~ dnorm(0,1e-2) }
}

initmod.pois.lv <- function() {
  ## Likelihood
  for(j in 1:s) { for(i in 1:n) { 
    eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
    mu[i,j] <- exp(eta[i,j])
    y[i,j] ~ dpois(mu[i,j])		
  } }
  
  for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } }
  
  ## Prior
  for(j in 1:s) { all.params[j,1] ~ dnorm(0,1e-2) } ## species-specific intercept	
  for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { 
    all.params[i,j] <- 0 } } ## Upper diagonal constrained to 0
  for(i in 1:num.lv) { all.params[i,i+1] ~ dunif(0,1e3) } ## Constraints on diagonal 	
  for(i in 2:num.lv) { for(j in 2:i) { 
    all.params[i,j] ~ dnorm(0,1e-2) } } ## Lower diagonal unconstrained
  for(i in (num.lv+1):s) { for(j in 2:(num.lv+1)) { 
    all.params[i,j] ~ dnorm(0,1e-2) } } ## All other elements
  for(j in 1:s) { all.params[j,num.lv+2] ~ dunif(0,1e3) } ## Remove this line for Poisson and Bernoulli
  
  for(i in 1:n) { site.params[i] ~ dnorm(0,1e-2) }
}

initmod.binom.lv <- function() {
  ## Likelihood
  for(j in 1:s) { for(i in 1:n) { 
    eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
    mu[i,j] <- exp(eta[i,j])/(1+exp(eta[i,j]))
    y[i,j] ~ dbern(mu[i,j])
  } }
  
  for(i in 1:n) { for(k in 1:num.lv) { lvs[i,k] ~ dnorm(0,1) } }
  
  ## Prior
  for(j in 1:s) { all.params[j,1] ~ dnorm(0,1e-2) } ## species-specific intercept	
  for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { all.params[i,j] <- 0 } } ## Upper diagonal constrained to 0
  for(i in 1:num.lv) { all.params[i,i+1] ~ dunif(0,1e3) } ## Constraints on diagonal 	
  for(i in 2:num.lv) { for(j in 2:i) { all.params[i,j] ~ dnorm(0,1e-2) } } ## Lower diagonal unconstrained
  for(i in (num.lv+1):s) { for(j in 2:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } ## All other elements
  for(j in 1:s) { all.params[j,num.lv+2] ~ dunif(0,1e3) } ## Remove this line for Poisson and Bernoulli
  
  for(i in 1:n) { site.params[i] ~ dnorm(0,1e-2) }
}

fit.initmod <- function(y, family, num.lv = 2, site.eff = FALSE, n.burnin = 5000, n.iter = 25000, n.thin = 5) {
  site.eff <- as.numeric(site.eff)
  n <- nrow(y); s <- ncol(y)
  if(num.lv >= ncol(y)) stop("# of LVs >= # of y...doesn't make sense.")
  
  jags.data <- list("y","n","s","num.lv","site.eff")
  jags.params <- c("all.params", "lvs")
  if(site.eff) { jags.params <- c(jags.params, "site.params") }
  
  if(family == "binomial") {
    jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = initmod.binom.lv, n.iter = n.iter, n.burnin = n.burnin, 
                    n.thin = n.thin, n.chains = 1, DIC = T, progress.bar = "text") }	
  if(family == "poisson") {
    jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = initmod.pois.lv, n.iter = n.iter, n.burnin = n.burnin, 
                    n.thin = n.thin, n.chains = 1, DIC = T, progress.bar = "text") }
  if(family == "negative.binomial") {
    jagsfit <- jags(data=jags.data, inits=NULL, jags.params, model.file = initmod.nb.lv, n.iter = n.iter, n.burnin = n.burnin, 
                    n.thin = n.thin, n.chains = 1, DIC = T, progress.bar = "text") }
  
  fit.mcmcBase <- jagsfit$BUGSoutput
  fit.mcmc <- mcmc(fit.mcmcBase$sims.matrix, start = 1, thin = 5) 
  
  if(family == "negative.binomial") { fit.mcmc[,(s*(num.lv+2)-s+1):(s*(num.lv+2))] <- 1/fit.mcmc[,(s*(num.lv+2)-s+1):(s*(num.lv+2))] }
  out.fit <- list(params.med = matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],2,median),nrow=s),
                  params.mean = matrix(apply(fit.mcmc[,grep("all.params",colnames(fit.mcmc))],2,mean),nrow=s), 
                  jags.model = jagsfit, site.eff = site.eff, family = family, 
                  lv.med = matrix(apply(fit.mcmc[,grep("lvs", colnames(fit.mcmc))],2,median),nrow=n),
                  lv.mean = matrix(apply(fit.mcmc[,grep("lvs", colnames(fit.mcmc))],2,mean),nrow=n))
  
  if(site.eff) {
    out.fit$site.params.med <- apply(fit.mcmc[,grep("site.params", colnames(fit.mcmc))],2,median)
    out.fit$site.params.mean <- apply(fit.mcmc[,grep("site.params", colnames(fit.mcmc))],2,mean) }
  if(!site.eff) { out.fit$site.params.med <- out.fit$site.params.mean <- rep(0,n) }
  
  out.fit$num.lv <- num.lv; 
  out.fit$y <- y; 
  out.fit$family <- family
  out.fit$site.eff <- site.eff
  out.fit$n.burnin <- n.burnin;
  out.fit$n.iter <- n.iter
  out.fit$n.thin <- n.thin
  
  return(out.fit) 
}	

############
## CORAL models for d > 1 and K > 1
############	
coralmod.nb.lv <- function() {
  ## Likelihood
  for(i in 1:N) { 
    clus.label[i] ~ dcat(mixprop[1:num.clus]) 
    for(k in 1:num.lv) { lvs[i,k] ~ dnorm(mu.clus[clus.label[i],k],1) } 	
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
      mu[i,j] <- exp(eta[i,j])
      u[i,j] ~ dgamma(1/all.params[j,num.lv+2],1/all.params[j,num.lv+2])
      y[i,j] ~ dpois(mu[i,j]*u[i,j]) ## Parameterizing the NB as a multiplicative random effect models
    } }
  
  ## Prior
  mixprop[1:num.clus] ~ ddirch(alpha[1:num.clus])
  for(k2 in 1:num.lv) { mu.clus[1,k2] <- -sum(mu.clus[2:num.clus,k2]) } ## Sum to zero constraint on cluster means
  for(k in 2:num.clus) { for(k2 in 1:num.lv) { mu.clus[k,k2] ~ dnorm(0,1e-2) } } 
  ## Identity covariance matrix
  for(k in 1:(num.lv-1)) { for(k2 in (k+1):num.lv) { Omega[k,k2] <- 0 } } 
  for(k in 2:num.lv) { for(k2 in 1:(k-1)) { Omega[k,k2] <- 0 } } 
  for(k in 1:num.lv) { Omega[k,k] <- 1 } 
  
  for(j in 1:S) { all.params[j,1] ~ dnorm(0,1e-2) } ## species-specific intercept	
  for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { all.params[i,j] <- 0 } } ## Upper diagonal constrained to 0
  for(i in 1:num.lv) { all.params[i,i+1] ~ dunif(0,1e3) } ## Constraints on diagonal
  for(i in 2:num.lv) { for(j in 2:i) { all.params[i,j] ~ dnorm(0,1e-2) } } ## Lower diagonal unconstrained
  for(i in (num.lv+1):S) { for(j in 2:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } ## All other elements
  for(j in 1:S) { all.params[j,num.lv+2] ~ dgamma(1,1) }
  
  for(i in 1:N) { site.params[i] ~ dnorm(0,1e-2) }
}

coralmod.pois.lv <- function() {
  ## Likelihood
  for(i in 1:N) { 
    clus.label[i] ~ dcat(mixprop[1:num.clus]) 
    for(k in 1:num.lv) { lvs[i,k] ~ dnorm(mu.clus[clus.label[i],k],1) } 	
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
      mu[i,j] <- exp(eta[i,j])
      y[i,j] ~ dpois(mu[i,j]) 
    } }
  
  ## Prior
  mixprop[1:num.clus] ~ ddirch(alpha[1:num.clus])
  for(k2 in 1:num.lv) { mu.clus[1,k2] <- -sum(mu.clus[2:num.clus,k2]) } ## Sum to zero constraint on cluster means
  for(k in 2:num.clus) { for(k2 in 1:num.lv) { mu.clus[k,k2] ~ dnorm(0,1e-2) } }
  ## Identity covariance matrix. 
  for(k in 1:(num.lv-1)) { for(k2 in (k+1):num.lv) { Omega[k,k2] <- 0 } } 
  for(k in 2:num.lv) { for(k2 in 1:(k-1)) { Omega[k,k2] <- 0 } } 
  for(k in 1:num.lv) { Omega[k,k] <- 1 } 
  
  for(j in 1:S) { all.params[j,1] ~ dnorm(0,1e-2) } ## species-specific intercept	
  for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { all.params[i,j] <- 0 } } ## Upper diagonal constrained to 0
  for(i in 1:num.lv) { all.params[i,i+1] ~ dunif(0,1e3) } ## constraints on diagonal
  for(i in 2:num.lv) { for(j in 2:i) { all.params[i,j] ~ dnorm(0,1e-2) } } ## Lower diagonal unconstrained
  for(i in (num.lv+1):S) { for(j in 2:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } ## All other elements
  
  for(i in 1:N) { site.params[i] ~ dnorm(0,1e-2) }
}

coralmod.binom.lv <- function() {
  ## Likelihood
  for(i in 1:N) { 
    clus.label[i] ~ dcat(mixprop[1:num.clus]) 
    for(k in 1:num.lv) { lvs[i,k] ~ dnorm(mu.clus[clus.label[i],k],1) } 	
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
      mu[i,j] <- exp(eta[i,j])/(1+exp(eta[i,j]))
      y[i,j] ~ dbern(mu[i,j]) 
    } }
  
  ## Prior
  mixprop[1:num.clus] ~ ddirch(alpha[1:num.clus])
  for(k2 in 1:num.lv) { mu.clus[1,k2] <- -sum(mu.clus[2:num.clus,k2]) } ## Sum to zero constraint on cluster means
  for(k in 2:num.clus) { for(k2 in 1:num.lv) { mu.clus[k,k2] ~ dnorm(0,1e-2) } }
  ## Identity covariance matrix. 
  for(k in 1:(num.lv-1)) { for(k2 in (k+1):num.lv) { Omega[k,k2] <- 0 } } 
  for(k in 2:num.lv) { for(k2 in 1:(k-1)) { Omega[k,k2] <- 0 } } 
  for(k in 1:num.lv) { Omega[k,k] <- 1 } 
  
  for(j in 1:S) { all.params[j,1] ~ dnorm(0,1e-2) } ## species-specific intercept	
  for(i in 1:(num.lv-1)) { for(j in (i+2):(num.lv+1)) { all.params[i,j] <- 0 } } ## Upper diagonal constrained to 0
  for(i in 1:num.lv) { all.params[i,i+1] ~ dunif(0,1e3) } ## Constraints on diagonal
  for(i in 2:num.lv) { for(j in 2:i) { all.params[i,j] ~ dnorm(0,1e-2) } } ## Lower diagonal unconstrained
  for(i in (num.lv+1):S) { for(j in 2:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } ## All other elements
  
  for(i in 1:N) { site.params[i] ~ dnorm(0,1e-2) }
}

############
## CORAL models for d = 1 and K > 1
############	
coralmod.nb.lv1 <- function() {
  ## Likelihood
  for(i in 1:N) { 
    clus.label[i] ~ dcat(mixprop[1:num.clus]) 
    for(k in 1:num.lv) { lvs[i,k] ~ dnorm(mu.clus[clus.label[i],k],1) }		
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
      mu[i,j] <- exp(eta[i,j])
      u[i,j] ~ dgamma(1/all.params[j,num.lv+2],1/all.params[j,num.lv+2])
      y[i,j] ~ dpois(mu[i,j]*u[i,j]) ## Parameterizing the NB as a multiplicative random effect models
    } }
  
  ## Prior
  mixprop[1:num.clus] ~ ddirch(alpha[1:num.clus])
  for(k2 in 1:num.lv) { mu.clus[1,k2] <- -sum(mu.clus[2:num.clus,k2]) } ## Sum to zero constraint on cluster means
  for(k in 2:num.clus) { for(k2 in 1:num.lv) { mu.clus[k,k2] ~ dnorm(0,1e-2) } }
  Omega <- 1	
  
  for(i in 1:S) { for(j in 1:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } 
  for(j in 1:S) { all.params[j,num.lv+2] ~ dgamma(0.01,0.01) }
  for(i in 1:N) { site.params[i] ~ dnorm(0,1e-2) }
}

coralmod.pois.lv1 <- function() {
  ## Likelihood
  for(i in 1:N) { 
    clus.label[i] ~ dcat(mixprop[1:num.clus]) 
    for(k in 1:num.lv) { lvs[i,k] ~ dnorm(mu.clus[clus.label[i],k],1) }		
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
      mu[i,j] <- exp(eta[i,j])
      y[i,j] ~ dpois(mu[i,j]) 
    } }
  
  ## Prior
  mixprop[1:num.clus] ~ ddirch(alpha[1:num.clus])
  for(k2 in 1:num.lv) { mu.clus[1,k2] <- -sum(mu.clus[2:num.clus,k2]) } ## Sum to zero constraint on cluster means
  for(k in 2:num.clus) { for(k2 in 1:num.lv) { mu.clus[k,k2] ~ dnorm(0,1e-2) } }
  Omega <- 1
  
  for(i in 1:S) { for(j in 1:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } 
  for(i in 1:N) { site.params[i] ~ dunif(-1e2,1e2) }
}

coralmod.binom.lv1 <- function() {
  ## Likelihood
  for(i in 1:N) { 
    clus.label[i] ~ dcat(mixprop[1:num.clus]) 
    for(k in 1:num.lv) { lvs[i,k] ~ dnorm(mu.clus[clus.label[i],k],1) }		
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i] + inprod(all.params[j,2:(num.lv+1)],lvs[i,])
      mu[i,j] <- exp(eta[i,j])/(1+exp(eta[i,j]))
      y[i,j] ~ dbern(mu[i,j]) 
    } }
  
  ## Prior
  mixprop[1:num.clus] ~ ddirch(alpha[1:num.clus])
  for(k2 in 1:num.lv) { mu.clus[1,k2] <- -sum(mu.clus[2:num.clus,k2]) } ## Sum to zero constraint on cluster means
  for(k in 2:num.clus) { for(k2 in 1:num.lv) { mu.clus[k,k2] ~ dnorm(0,1e-2) } }
  Omega <- 1
  
  for(i in 1:S) { for(j in 1:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } 
  for(i in 1:N) { site.params[i] ~ dunif(-1e2,1e2) }
}

############
## CORAL models for d = 0 and K = 0
############     
coralmod.nb.lv0clus0 <- function() {
  ## Likelihood
  for(i in 1:N) { 
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i]
      mu[i,j] <- exp(eta[i,j])
      u[i,j] ~ dgamma(1/all.params[j,2],1/all.params[j,2])
      y[i,j] ~ dpois(mu[i,j]*u[i,j]) ## Parameterizing the NB as a multiplicative random effect models
    } }
  
  ## Prior
  Omega <- 1
  for(i in 1:S) { for(j in 1:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } 
  for(i in 1:N) { site.params[i] ~ dnorm(0,1e-2) }
  for(j in 1:S) { all.params[j,2] ~ dgamma(0.1,0.1) } ## Stronger prior needed to prevent problems with sample in null models
}

coralmod.pois.lv0clus0 <- function() {
  ## Likelihood
  for(i in 1:N) { 
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i]
      mu[i,j] <- exp(eta[i,j])
      y[i,j] ~ dpois(mu[i,j]) 
    } }
  
  ## Prior
  Omega <- 1     
  for(i in 1:S) { for(j in 1:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } 
  for(i in 1:N) { site.params[i] ~ dunif(-1e2,1e2) }
}

coralmod.binom.lv0clus0 <- function() {
  ## Likelihood
  for(i in 1:N) { 
    for(j in 1:S) {
      eta[i,j] <- all.params[j,1] + site.eff*site.params[i]
      mu[i,j] <- exp(eta[i,j])/(1+exp(eta[i,j]))
      y[i,j] ~ dbern(mu[i,j]) 
    } }
  
  ## Prior
  Omega <- 1     
  for(i in 1:S) { for(j in 1:(num.lv+1)) { all.params[i,j] ~ dnorm(0,1e-2) } } 
  for(i in 1:N) { site.params[i] ~ dnorm(0,1e-2) }
}

############
## Auxilary functions
#############
## Density for a mixture of multivariate gaussians
dlvmixture <- function(b.i, pis, mu.clus, Omega) {
  lik <- 0; K <- length(pis); 
  for(k in 1:K) { lik <- lik + pis[k]*dmvnorm(b.i,mu.clus[k,],Omega) }
  lik }

## Obtain the mode in each row of matrix; used to find the group-label with highest posterior classification probability
get.mode <- function(x) {
  ux <- unique(x); 
  ux[which.max(tabulate(match(x, ux)))] }

## Calculate posterior classification probabilities
calc.taus <- function(data, mixprops, mu.clus, lvs) {
  taus <- calc.f <- matrix(NA,nrow(data),length(mixprops))
  for(i in 1:nrow(data)) { 
    for(k in 1:length(mixprops)) { calc.f[i,k] <- taus[i,k] <- mixprops[k]*dmvnorm(lvs[i,], mu.clus[k,], diag(rep(1,length(lvs[i,])))) }
    taus[i,] <- taus[i,]/sum(taus[i,]) }
  return(list(taus = taus, calc.f = calc.f)) }

## Calculate conditional log-likelihood: loglik = sum_{i=1}^n \sum_{j=1}^s log f(y_ij|b_i) 
get.cond.predlogl <- function(y, lvs, params, site.params = NULL, family, phi = T) {
  n <- nrow(y); s <- ncol(y); 
  num.lv <- ncol(as.matrix(lvs)); 
  lvs <- as.matrix(lvs)
  
  predlogl.mat <- fitted.vals <- matrix(NA,n,s) 
  if(is.null(site.params)) site.params <- rep(0,n)
  
  for(i in 1:n) { 
    eta <- params[,1:(num.lv+1)]%*%c(1,lvs[i,]) + site.params[i]; #print(eta)
    if(family == "negative.binomial") {
      fitted.vals[i,] <- exp(eta)
      if(phi) { predlogl.mat[i,] <- dnbinom(as.vector(unlist(y[i,])), mu = exp(eta), size = 1/params[,ncol(params)], log = T) }
      if(!phi) { predlogl.mat[i,] <- dnbinom(as.vector(unlist(y[i,])), mu = exp(eta), size = params[,ncol(params)], log = T) }
    }
    if(family == "binomial") {
      fitted.vals[i,] <- exp(eta)/(1+exp(eta))
      predlogl.mat[i,] <- dbinom(as.vector(unlist(y[i,])), 1, p = exp(eta)/(1+exp(eta)), log = T) }
    if(family == "poisson") {
      fitted.vals[i,] <- exp(eta)
      predlogl.mat[i,] <- dpois(as.vector(unlist(y[i,])), lambda = exp(eta), log = T) }
  }
  
  return(list(predlogl.mat = predlogl.mat, fitted.vals = fitted.vals)) 
}

## Calculate marginal likelihood using importance sampling: loglik = sum_{i=1}^n \log( \int \prod_{j=1}^s f(y_ij|b_i) f(b_i) db_i )
get.marg.predlogl <- function(y, params, site.params, pis, mu.clus, family) {
  n <- nrow(y); s <- ncol(y); 
  num.lv <- ncol(mu.clus); 
  num.clus <- length(pis)
  loglik <- 0; reps <- 20000; 
  Omega <- diag(num.lv)
  
  big.gen.lv <- rmvnorm(reps,mean=colMeans(mu.clus)) ## Use multivariate normal with mean vector each to average of mu.clus as importance density -- crude but OK
  dmix <- dimpor <- numeric(reps)
  dimpor <- dmvnorm(big.gen.lv,mean=colMeans(mu.clus)) 
  for(t in 1:reps) dmix[t] <- dlvmixture(big.gen.lv[t,], pis, mu.clus, Omega)
  
  site.logl.reps <- matrix(0,reps,n)
  for(i in 1:n) { for(t in 1:reps) {
    site.logl.reps[t,i] <- exp(apply(get.cond.predlogl(matrix(y[i,],1), lvs = big.gen.lv[t,], params, site.params = site.params[i], family, phi = T)$predlogl.mat,1,sum)) 
  } }
  site.logl.reps[!is.finite(site.logl.reps)] <- 1
  
  out <- log(apply(site.logl.reps*dmix/dimpor,2,mean))
  out }

## Fitted values and linear predictors
fitted.coral <- function(coralfit, est = "median") {
  y <- coralfit$y
  family <- coralfit$family
  n <- nrow(y); s <- ncol(y); 
  num.lv <- ncol(coralfit$lv.med) 
  fitted.out <- matrix(NA,n,s)
  
  if(is.null(coralfit$lv.med)) { eta <- matrix(1,n,1)%*%t(coralfit$params.med[,1:(num.lv+1)]) }
  if(!is.null(coralfit$lv.med)) { eta <- cbind(1,coralfit$lv.med)%*%t(coralfit$params.med[,1:(num.lv+1)]) }
  if(est == "mean") {
    if(is.null(coralfit$lv.mean)) { eta <- matrix(1,n,1)%*%t(coralfit$params.mean[,1:(num.lv+1)]) }
    if(!is.null(coralfit$lv.mean)) { eta <- cbind(1,coralfit$lv.mean)%*%t(coralfit$params.mean[,1:(num.lv+1)]) 
    } }
  
  for(i in 1:n) {
    if(est == "median") { if(!is.null(coralfit$site.params.med)) eta[i,] <- eta[i,] + coralfit$site.params.med[i] }
    if(est == "mean") { if(!is.null(coralfit$site.params.mean)) eta[i,] <- eta[i,] + coralfit$site.params.mean[i] 
    } }
  
  for(i in 1:n) {
    if(family == "binomial") fitted.out[i,] <- exp(eta[i,])/(1+exp(eta[i,]))
    if(family == "poisson") fitted.out[i,] <- exp(eta[i,])
    if(family == "negative.binomial") fitted.out[i,] <- exp(eta[i,]) }
  
  return(list(fitted = fitted.out, lin.pred = eta))
}

## Dunn-Smyth residuals
ds.residuals <- function(coralfit, est = "median") {  
  y <- coralfit$y
  family <- coralfit$family
  n <- nrow(y); s <- ncol(y); 
  num.lv <- ncol(coralfit$lv.med) 
  mus <- fitted.coral(y, coralfit, family, est)$fitted
  
  if(est == "median" & family == "negative.binomial") phis <- coralfit$params.median[,num.lv+2]
  if(est == "mean" & family == "negative.binomial") phis <- coralfit$params.mean[,num.lv+2]
  
  ds.res.out <- y
  for(i in 1:n) { for(j in 1:s) {
    if(family == "poisson") { a <- ppois(as.vector(unlist(y[i,j]))-1, mus[i,j]); b <- ppois(as.vector(unlist(y[i,j])), mus[i,j]) }
    if(family == "negative.binomial") { a <- pnbinom(as.vector(unlist(y[i,j]))-1, mu=mus[i,j], size=1/phis[j]); b <- pnbinom(as.vector(unlist(y[i,j])), mu=mus[i,j], size=1/phis[j]) }
    if(family == "binomial") { a <- pbinom(as.vector(unlist(y[i,j]))-1, 1, mus[i,j]); b <- pbinom(as.vector(unlist(y[i,j])), 1, mus[i,j]) }
    
    u <- runif(n = 1, min = a, max = b)
    ds.res.out[i,j] <- qnorm(u) 
  } }
  
  ds.res.out
}
