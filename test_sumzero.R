test_code <- nimbleCode({
  for(coef in 2:p){
    beta[coef] ~ dnorm(0, sd = 10)
  }
  beta[1] <- -sum(beta[2:p])
  beta0 ~ dnorm(0, sd = 10)
  sigma ~ dunif(0, 100)
  
  for(i in 1:n){
    y[i] ~ dnorm(beta0 + inprod(beta[1:p], x[i, 1:p]), sd = sigma)
  }
})

# generate data
set.seed(1)
n <- 1000
p <- 5
X <- matrix(rnorm(n*p), n, p)
beta0 <- rnorm(1)
beta <- rnorm(p)
beta[1] <- -sum(beta[2:p])
y <- rnorm(n, beta0 + c(X%*%beta), sd = 1)

beta_init <- rnorm(5)
beta_init[1] <- -sum(beta_init[2:5])

test <- fit_model(
  seed = 1,
  code = test_code,
  data = list(
    y = y,
    x = X
  ),
  constants = list(
    n = n,
    p = p
  ),
  inits = list(
    beta0 = rnorm(1),
    beta = beta_init, 
    sigma = abs(rnorm(1))
  ),
  niter = 10000,
  burnin = 5000,
  nchains = 3,
  thin = 1
)
colMeans(
  do.call("rbind", test)
)
c(beta, beta0, 1)

sum(test[[1]][2,1:5])








