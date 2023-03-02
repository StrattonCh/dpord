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
  # for(dim in 1:d){
  #   mu[1, dim] <- -sum(mu[2:K, dim])
  # }
  for(clus in 1:K){
    for(dim in 1:d){
      mu[clus, dim] ~ dnorm(0, var = 10)
    }
  }
  # for(i in 1:K){
  #   mu[i, 1:d] ~ dmnorm(mu0[1:d], Lambda0[1:d, 1:d])
  # }
  
  
  for(site in 1:nsites){
    for(dim in 1:d){
      z[site, dim] ~ dnorm(mu[clus_id[site], dim], var = 1)
    }
    # # identity matrix for constraint
    # z[site, 1:d] ~ dmnorm(mu[clus_id[site], 1:d], cov = S[1:d, 1:d])
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
init <- init_func(mat, d = 2, max_clus = 2)
init$clus_id <- rep(1, nrow(mat))
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
  niter = niter,
  burnin = nburnin,
  nchains = 1,
  thin = thin
)
