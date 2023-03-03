y <- c(
  rnorm(30, -3, .5),
  rnorm(30, 0, .5),
  rnorm(30, 3, .5)
)
hist(y)

test_code <- nimbleCode({
  for(i in 1:K){
    mu[i] ~ dnorm(0, sd = 10)
  }
  # mu[1] <- -sum(mu[2:K])
  
  for(i in 1:n){
    clus_id[i] ~ dcat(pi[1:K])
  }
  pi[1:K] ~ ddirch(ones[1:K])
  sigma ~ dunif(0, 20)
  
  for(i in 1:n){
    y[i] ~ dnorm(mu[clus_id[i]], sigma)
  }
})

inits <- list(
  mu = rnorm(3),
  clus_id = sample(1:3, size = length(y), replace = T),
  pi = rep(1/3, 3),
  sigma = 1
)
inits$mu[1] <- -sum(inits$mu[2:3])
model <- nimbleModel(
  code = test_code,
  data = list(
    y = y
  ),
  constants = list(
    n = length(y),
    K = 3,
    ones = rep(1, 3)
  ),
  init = list(
    mu = rnorm(3)
  )
)

cmodel <- compileNimble(model)
rMCMC <- buildMCMC(cmodel)
cMCMC <- compileNimble(rMCMC)
results <- runMCMC(cMCMC, niter = 10000, nburnin = 5000)
