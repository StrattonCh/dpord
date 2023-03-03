coral_code <- nimbleCode({
  # priors
  phi ~ dunif(0, 10)
  ## site effects
  for(site in 1:nsites){
    alpha[site] ~ dnorm(0, 1)
  }
  
  ## species effects
  for(species in 1:nspecies){
    beta[species] ~ dnorm(0, 1)
  }
  
  ## z priors
  ### finite mixture
  for(i in 1:nsites){
    clus_id[i] ~ dcat(pi[1:K])
  }
  pi[1:K] ~ ddirch(ones[1:K])
  
  ### table parameters - fix covariance as identity
  mu1[1] <- -sum(mu1[2:K])
  mu2[1] <- -sum(mu2[2:K])
  for(clus in 2:K){
    mu1[clus] ~ dnorm(0, var = 10)
    mu2[clus] ~ dnorm(0, var = 10)
  }

  
  
  for(site in 1:nsites){
    z[site,1] ~ dnorm(mu1[clus_id[site]], 1)
    z[site,2] ~ dnorm(mu2[clus_id[site]], 1)
  }
  
  # theta prior
  ## upper triangle = 0
  for(row in 1:(d-1)){
    for(col in (row+1):d){
      theta[row, col] <- 0
    }
  }
  
  ## diag > 0
  for(diag_element in 1:d){
    theta[diag_element, diag_element] ~ T(dnorm(0, sd = 1), 0, Inf)
  }
  
  ## lower diag of first d rows
  for(row in 2:d){
    for(col in 1:(row-1)){
      theta[row, col] ~ dnorm(0, sd = 1)
    }
  }
  
  ## all other elements
  for(row in (d+1):nspecies){
    for(col in 1:d){
      theta[row, col] ~ dnorm(0, sd = 1)
    }
  }
  
  # likelihood 
  for(site in 1:nsites){
    for(species in 1:nspecies){
      log(lambda[site, species]) <- alpha[site] + beta[species] + inprod(z[site,1:d], theta[species, 1:d])
      nbp[site, species] <- phi / (phi + lambda[site, species])
      Y[site, species] ~ dnegbin(prob = nbp[site, species], size = phi)
    }
  }
})
init_func <- function(mat, d = 2, max_clus){
  nrow <- nrow(mat)
  nspecies <- ncol(mat)
  theta_init <- matrix(rnorm(d * nspecies), nspecies, d)
  diag(theta_init) <- abs(rnorm(d))
  theta_init[upper.tri(theta_init)] <- 0
  mu <- matrix(rnorm(d*max_clus), max_clus, d)
  if(max_clus != 1) mu[1,] <- -colSums(mu[-1,,drop = F])
  
  list(
    alpha = rnorm(nrow),
    beta = rnorm(nspecies),
    z = matrix(rnorm(nrow * d), nrow, d),
    theta = theta_init,
    dp_con = 1,
    clus_id = sample(c(1, 2), size = nrow, replace = T),
    mu = mu,
    phi = runif(1)
  )
}

# fit model
ndx <- 2

init <- list(
  alpha = rnorm(30),
  beta = rnorm(27),
  mu1 = rnorm(2),
  mu2 = rnorm(2),
  clus_id = sample(c(1:2), 30, replace = T)
)


test <- fit_model(
  seed = 1,
  code = coral_code,
  data = list(
    Y = mat
  ),
  constants = list(
    nsites = nrow(mat),
    nspecies = ncol(mat),
    d = 2,
    K = ndx,
    S = diag(2),
    mu0 = rep(0, 2),
    Lambda0 = diag(2),
    ones = rep(1, ndx)
  ),
  inits = init,
  niter = 50000,
  burnin = 25000,
  nchains = 1,
  thin = 5
)
colnames(test)
test[,85]

# step by step
model <- nimbleModel(
  code = coral_code,
  data = list(
    Y = mat
  ),
  constants = list(
    nsites = nrow(mat),
    nspecies = ncol(mat),
    d = 2,
    K = ndx,
    S = diag(2),
    mu0 = rep(0, 2),
    Lambda0 = diag(2),
    ones = rep(1, ndx)
  ),
  inits = init
)

cmodel <- compileNimble(model)


