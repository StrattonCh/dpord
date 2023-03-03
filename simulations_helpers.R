# helper functions for plotting, etc
plot_true_and_realizations <- function(plot_grid, compas_params, ndx = 1){
  # function to panel true mixture distribution with realizations
  
  # housekeeping
  params <- compas_params[[ndx]]
  k <- length(params$sigma)
  tmp <- sapply(1:k, function(x){
    mvtnorm::dmvnorm(plot_grid, mean = params$mu[x,], sigma = params$sigma[x]*diag(2))
  })
  
  # true density
  p1 <- tibble(
    prob = c(tmp),
    group = factor(rep(1:k, each = nrow(plot_grid))),
    x1 = rep(plot_grid[,1], k),
    x2 = rep(plot_grid[,2], k)
  ) %>%
    ggplot() + 
    geom_contour(
      aes(x = x1, y = x2, z = prob, col = group)
    ) +
    theme_bw() +
    labs(title = paste0("True mixture, k = ", k))
  
  # realizations
  out <- list()
  for(i in 1:8){
    tmp <- lapply(1:k, function(x){
      mvtnorm::rmvnorm(15, mean = params$mu[x,], sigma = params$sigma[x]*diag(2))
    })
    tmp <- do.call("rbind", tmp)
    
    out[[i]] <- tibble(
      group = factor(rep(1:k, each = 15)),
      x1 = tmp[,1],
      x2 = tmp[,2]
    ) %>%
      ggplot() +
      geom_point(aes(x = x1, y = x2, col = group)) +
      theme_bw() +
      labs("Realization of true mixture")
  }
  
  gridExtra::grid.arrange(
    p1, out[[1]],out[[2]],out[[3]],
    out[[4]],out[[5]],out[[6]],
    out[[7]],out[[8]], nrow = 3 
  )
  
}
nimble_summary <- function(fit, warmup = nrow(fit[[1]])/2, thin = 1){
  # convert to coda for normal summary
  fit_warmup <- lapply(fit, function(x) x[(warmup+1):nrow(x),])
  coda_samples <- as.mcmc.list(lapply(fit_warmup, function(x) as.mcmc(
    x, start = warmup+1, end = nrow(fit), thin = thin
  )))
  
  sum <- summary(coda_samples)
  params <- dimnames(sum$statistics)[[1]]
  tmp_sum <- cbind(sum$statistics, sum$quantiles)
  
  # get r hat / n_eff
  mat <- matrix(NA, nrow = nrow(tmp_sum), ncol = 3)
  colnames(mat) <- c("Rhat", "ess_bulk", "ess_tail")
  for(i in 1:nrow(tmp_sum)){
    tmp <- sapply(fit, function(x) x[,i])
    mat[i,] <- c(Rhat(tmp), ess_bulk(tmp), ess_tail(tmp))
  }
  
  # out 
  out <- cbind(tmp_sum, mat)
  return(out)
}
waic_dpord <- function(Y, mcmc, frame_ind = NULL, hier = F){
  if(hier){
    # get params 
    mu <- mcmc$mu
    clus_id <- mcmc$clus_id
    other <- mcmc$other
    alpha <- other[,grepl("alpha", colnames(other))]
    beta <- other[,grepl("beta", colnames(other))]
    gamma <- other[,grepl("gamma[[]", colnames(other))]
    z <- other[,grepl("z", colnames(other))]
    theta <- other[,grepl("theta", colnames(other))]
    
    # get dimension
    d <- sapply(colnames(mu), function(x){
      str_sub(x, start = nchar(x)-1, end = nchar(x)-1)
    }) %>% as.numeric %>% max
    
    # compute lppd
    ## lik = alpha[i] + beta[j] + z[i] theta[j] + gamma[i]
    log_py <- matrix(NA, nrow = nrow(z), ncol = length(c(Y)))
    message("Calculating log likelihood")
    pb <- txtProgressBar(min = 0, max = nrow(log_py), style = 3, char = "=")
    for(iter in 1:nrow(log_py)){
      alpha_vec_s <- alpha[iter,]
      beta_vec_s <- beta[iter,]
      z_matrix_s <- matrix(z[iter,], ncol = d)
      z_matrix_s <- z_matrix_s[rep(1:nrow(z_matrix_s), table(frame_ind)),]
      theta_matrix_s <- matrix(theta[iter,], ncol = d)
      gamma_s <- gamma[iter, ]
      
      linpred <- rep(rep(alpha_vec_s, table(frame_ind)), ncol(Y)) +
        rep(beta_vec_s, each = nrow(Y)) +
        c(z_matrix_s %*% t(theta_matrix_s)) +
        rep(gamma_s, ncol(Y))
      pi <- unname(exp(linpred) / (1 + exp(linpred)))
      
      # fix numerical issues
      pi[which(pi > 0.99999)] <- 0.99999
      pi[which(pi < 0.00001)] <- 0.00001
      
      log_py[iter,] <- c(Y) * log(pi) + (1 - c(Y)) * log(1 - pi)
      setTxtProgressBar(pb, iter)
    }
    close(pb)
    
    # compute lppd
    message("Computing lppd")
    lppd <- sum(log(colMeans(exp(log_py))))
    
    # compute penalty
    message("Computing penalty")
    pwaic2 <- sum(apply(log_py, 2, function(x) var(x)))
    
    # compute waic
    waic <- -2*(lppd - pwaic2)
    out <- waic
  } else{
    # get params 
    mu <- mcmc$mu
    clus_id <- mcmc$clus_id
    other <- mcmc$other
    alpha <- other[,grepl("alpha", colnames(other))]
    beta <- other[,grepl("beta", colnames(other))]
    z <- other[,grepl("z", colnames(other))]
    theta <- other[,grepl("theta", colnames(other))]
    
    # get dimension
    d <- sapply(colnames(mu), function(x){
      str_sub(x, start = nchar(x)-1, end = nchar(x)-1)
    }) %>% as.numeric %>% max
    
    # compute lppd
    ## lik = alpha[i] + beta[j] + z[i] theta[j] + gamma[i]
    log_py <- matrix(NA, nrow = nrow(z), ncol = length(c(Y)))
    message("Calculating log likelihood")
    for(iter in 1:nrow(log_py)){
      alpha_vec_s <- alpha[iter,]
      beta_vec_s <- beta[iter,]
      z_matrix_s <- matrix(z[iter,], ncol = d)
      theta_matrix_s <- matrix(theta[iter,], ncol = d)
      
      linpred <- rep(alpha_vec_s, ncol(Y)) +
        rep(beta_vec_s, each = nrow(Y)) +
        c(z_matrix_s %*% t(theta_matrix_s)) 
      pi <- unname(exp(linpred) / (1 + exp(linpred)))
      
      # fix numerical issues
      pi[which(pi > 0.99999)] <- 0.99999
      pi[which(pi < 0.00001)] <- 0.00001
      
      log_py[iter,] <- c(Y) * log(pi) + (1 - c(Y)) * log(1 - pi)
    }
    
    # compute lppd
    message("Computing lppd")
    lppd <- sum(log(colMeans(exp(log_py))))
    
    # compute penalty
    message("Computing penalty")
    pwaic2 <- sum(apply(log_py, 2, function(x) var(x)))
    
    # compute waic
    waic <- -2*(lppd - pwaic2)
    out <- waic
  }
  
  return(out)
}
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# dpord functions
dpord_code <- nimbleCode({
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
  ### Dirichlet process mixture parameters
  clus_id[1:nsites] ~ dCRP(dp_con, size = nsites)
  dp_con ~ dgamma(1, 2)
  
  ### table parameters - fix covariance as identity
  for(i in 1:max_clus){
    mu[i, 1:d] ~ dmnorm(mu0[1:d], Lambda0[1:d, 1:d])
  }
  
  for(site in 1:nsites){
    # identity matrix for constraint
    z[site, 1:d] ~ dmnorm(mu[clus_id[site], 1:d], cov = S[1:d, 1:d])
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
  for(i in 1:K){
    mu[i, 1:d] ~ dmnorm(mu0[1:d], Lambda0[1:d, 1:d])
  }
  
  for(site in 1:nsites){
    # identity matrix for constraint
    z[site, 1:d] ~ dmnorm(mu[clus_id[site], 1:d], cov = S[1:d, 1:d])
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
coral_code_k1 <- nimbleCode({
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
  
  for(site in 1:nsites){
    # identity matrix for constraint
    z[site, 1:d] ~ dmnorm(mu0[1:d], cov = S[1:d, 1:d])
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
  max_clus <- nrow
  theta_init <- matrix(rnorm(d * nspecies), nspecies, d)
  diag(theta_init) <- abs(rnorm(d))
  theta_init[upper.tri(theta_init)] <- 0
  mu <- matrix(rnorm(d*max_clus), max_clus, d)
  if(max_clus != 1) mu[1,] <- -colSums(mu[-1,])
  
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
fit_model <- function(seed = 1, code, data, constants, inits, niter, nchains, thin = 1, burnin = 0, addMoni = c("clus_id", "z")){
  library(nimble)
  
  # R model
  model <- nimbleModel(code, constants, data)
  
  # C model
  model_c <- compileNimble(model)
  
  # R mcmc
  model_conf <- configureMCMC(model)
  if(!is.null(addMoni)){
    model_conf$addMonitors(addMoni)
  }
  
  
  # R mcmc
  mcmc <- buildMCMC(model_conf)
  
  # C mcmc
  mcmc_c <- compileNimble(mcmc, project = model_c)
  
  # run model
  out <- runMCMC(
    mcmc_c, 
    niter = niter, 
    nchains = nchains, 
    thin = thin, 
    init = inits,
    nburnin = burnin,
    setSeed = seed
  )
  
  # out
  return(out)
}
fit_model_k1 <- function(seed = 1, code, data, constants, inits, niter, nchains, thin = 1, burnin = 0){
  library(nimble)
  
  # R model
  model <- nimbleModel(code, constants, data)
  
  # C model
  model_c <- compileNimble(model)
  
  # R mcmc
  model_conf <- configureMCMC(model)
  model_conf$addMonitors(c("z"))
  
  # R mcmc
  mcmc <- buildMCMC(model_conf)
  
  # C mcmc
  mcmc_c <- compileNimble(mcmc, project = model_c)
  
  # run model
  out <- runMCMC(
    mcmc_c, 
    niter = niter, 
    nchains = nchains, 
    thin = thin, 
    init = inits,
    nburnin = burnin,
    setSeed = seed
  )
  
  # out
  return(out)
}
ord_ls <- function(mcmc, d = 2, seed = NULL, force_K = NULL){
  # optional seed
  if(!is.null(seed)) set.seed(seed)
  
  # housekeeping
  mu <- lapply(mcmc, function(x) x[,which(grepl("mu", colnames(x)))])
  clus_id <- lapply(mcmc, function(x) x[,which(grepl("clus_id", colnames(x)))])
  other_ndx <- c(1:length(colnames(mcmc[[1]])))[-c(
    sort(c(
      which(grepl("clus_id", colnames(mcmc[[1]]))), 
      which(grepl("mu", colnames(mcmc[[1]])))
    ))
  )]
  other <- lapply(mcmc, function(x) x[,other_ndx])
  
  mu_mcmc <- do.call("rbind", mu)
  clus_mcmc <- do.call("rbind", clus_id)
  other_mcmc <- do.call("rbind", other)
  
  # helper function
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  # step 1 - determine number of non-empty cluster for each iteration
  n_clus_mcmc <- apply(clus_mcmc, 1, function(x) length(unique(x)))
  
  # step 2 - estimate the mode
  if(is.null(force_K)){
    n_clus_mode <- getmode(n_clus_mcmc)
  } else{
    n_clus_mode <- force_K
  }
  
  # step 3 - filter to iterations with n_clus_mcmc = n_clus_mode
  ## filter first
  which_equal_mode_index <- which(n_clus_mcmc == n_clus_mode)
  mu_mcmc_filtered <- mu_mcmc[which_equal_mode_index, ]
  clus_mcmc_filtered <- clus_mcmc[which_equal_mode_index, ]
  other_mcmc_filtered <- other_mcmc[which_equal_mode_index, ]
  
  ## remove empty clusters - loop
  mu_reduced <- matrix(NA, nrow = nrow(mu_mcmc_filtered), ncol = n_clus_mode * d)
  clus_reduced <- array(NA, dim = dim(clus_mcmc_filtered))
  for(i in 1:length(which_equal_mode_index)){
    clus_id_i <- sort(unique(clus_mcmc_filtered[i, ]))
    col_ndx <- which(
      stringr::str_sub(
        colnames(mu_mcmc_filtered), 
        start = stringr::str_locate(colnames(mu_mcmc_filtered), pattern = "\\[")[,1] + 1,
        end = stringr::str_locate(colnames(mu_mcmc_filtered), pattern = "\\,")[,1] - 1
      ) %in% clus_id_i
    )
    
    mu_reduced[i, ] <- mu_mcmc_filtered[i,col_ndx]
    colnames(mu_reduced) <- paste0("mu[", rep(1:n_clus_mode, d), ", ", rep(1:d, each = n_clus_mode), "]")
    clus_reduced[i, ] <- as.numeric(factor(clus_mcmc_filtered[i,]))
    colnames(clus_reduced) <- colnames(clus_mcmc_filtered)
  }
  
  # step 4 - create data matrix, cluster with k-means
  data_matrix <- matrix(
    c(mu_reduced), 
    nrow = nrow(mu_reduced) * n_clus_mode,
    ncol = d
  )
  k_means <- kmeans(data_matrix, centers = n_clus_mode)
  rho <- matrix(k_means$cluster, nrow = nrow(mu_reduced), ncol = n_clus_mode)
  
  # step 5 - check whether rho_i is a permutation of 1:n_clus_mode
  keep_ndx <- which(apply(rho, 1, function(x) all(sort(x) == 1:n_clus_mode)))
  mu_reduced <- mu_reduced[keep_ndx,]
  clus_reduced <- clus_reduced[keep_ndx,]
  other_reduced <- other_mcmc_filtered[keep_ndx, ]
  rho <- rho[keep_ndx,] %>% as.matrix(., ncol = n_clus_mode)
  
  # step 6 - relabel according to rho
  ## relabel mus first - works, but slow
  mu_relabeled <- array(NA, dim = dim(mu_reduced))
  colnames(mu_relabeled) <- colnames(mu_reduced)
  for(i in 1:nrow(mu_relabeled)){
    for(j in 1:ncol(mu_relabeled)){
      # get group id
      group_id <- stringr::str_sub(
        colnames(mu_reduced)[j], 
        start = stringr::str_locate(colnames(mu_reduced)[j], pattern = "\\[")[,1] + 1,
        end = stringr::str_locate(colnames(mu_reduced)[j], pattern = "\\,")[,1] - 1
      ) %>% as.numeric()
      
      dim_id <- stringr::str_sub(
        colnames(mu_reduced)[j], 
        start = stringr::str_locate(colnames(mu_reduced)[j], pattern = "\\,")[,1] + 1,
        end = stringr::str_locate(colnames(mu_reduced)[j], pattern = "\\]")[,1] - 1
      ) %>% as.numeric()
      
      # replace with proper value according to rho
      rho_i <- rho[i, ]
      mu_relabeled[i,j] <- mu_reduced[i, which(
        colnames(mu_reduced) == paste0("mu[", which(rho_i == group_id), ", ", dim_id, "]")
      )]
    }
  }
  
  ## relabel clus_id's
  clus_relabeled <- array(NA, dim = dim(clus_reduced))
  colnames(clus_relabeled) <- colnames(clus_reduced)
  for(i in 1:nrow(clus_relabeled)){
    rho_i <- rho[i,]
    clus_id_i <- clus_reduced[i,]
    clus_relabeled[i,] <- sapply(clus_id_i, function(x) rho_i[x])
  }
  
  out <- list(
    mu = mu_relabeled, 
    clus_id = clus_relabeled, 
    other = other_reduced,
    which_iter = which_equal_mode_index[keep_ndx]
  )
  
  return(out)
}
summarize_dpord <- function(fit, seed = 1){
  # housekeeping
  function_out <- list()
  nchains <- length(fit)
  niter <- nrow(fit[[1]])
  function_out$niter <- niter
  function_out$nchains <- nchains
  
  # label-switching algorithm first
  samples_rl <- ord_ls(mcmc = fit, d = 2, seed = seed)
  
  # grab pseudo convergence diagnostics
  # these are meant for balanced chains
  # subset to only iterations present across all chains
  out <- list()
  for(chain in 1:nchains){
    out[[chain]] <- with(samples_rl, {
      which_iter[
        (which_iter <= niter*chain) &
          (which_iter >= niter*(chain-1))
      ]
    })
  }
  ndx <- sapply(out, function(x) x[1:min(sapply(out, length))])
  if(any(sapply(out, length) < 100)) return(NA)
  ndx <- ndx[-c(1, nrow(ndx)),]
  all_samples <- cbind(samples_rl$mu, samples_rl$other)
  all_samples <- all_samples[,-which(colnames(all_samples) == "theta[1, 2]")]
  conv_diag <- matrix(NA, nrow = ncol(all_samples), ncol = 3)
  for(col in 1:ncol(all_samples)){
    tmp <- all_samples[which(samples_rl$which_iter %in% c(ndx)), col]
    tmp2 <- matrix(tmp, nrow = nrow(ndx))
    conv_diag[col, 1:3] <- c(
      Rhat(tmp2), ess_bulk(tmp2), ess_tail(tmp2)
    )
  }
  conv_diag_tbl <- tibble(
    param = colnames(all_samples),
    rhat = conv_diag[,1],
    ess_bulk = conv_diag[,2],
    ess_tail = conv_diag[,3]
  )
  function_out$conv <- conv_diag_tbl
  
  # number of clusters based on raw mcmc
  clus_id <- lapply(fit, function(x) x[,which(grepl("clus_id", colnames(x)))])
  clus_id_all <- do.call("rbind", clus_id)
  ngroups_mcmc <- apply(clus_id_all, 1, function(x) length(unique(x)))
  tab <- round(table(ngroups_mcmc) / length(ngroups_mcmc), 3)
  function_out$raw_nclus <- as.numeric(names(which(tab == max(tab))))
  
  # modal cluster assignments and means for z's
  z_mcmc <- samples_rl$other[,which(grepl("z", colnames(samples_rl$other)))]
  z_mcmc_dim1 <- z_mcmc[,grepl(", 1", colnames(z_mcmc))]
  z_mcmc_dim2 <- z_mcmc[,grepl(", 2", colnames(z_mcmc))]
  z_tbl <- tibble(
    obs_id = 1:ncol(z_mcmc_dim1),
    z1 = colMeans(z_mcmc_dim1),
    z2 = colMeans(z_mcmc_dim2),
    modal_clus = apply(samples_rl$clus_id, 2, getmode),
    clus_prob = colMeans(
      samples_rl$clus_id == matrix(rep(apply(samples_rl$clus_id, 2, getmode), each = nrow(samples_rl$clus_id)), nrow = nrow(samples_rl$clus_id))
    )
  )
  function_out$z_tbl <- z_tbl
  function_out$modal_nclus <- length(unique(z_tbl$modal_clus))
  function_out$nrow_relabeled <- length(samples_rl$which_iter)

  # out
  return(function_out)
}
summarize_coral <- function(fit, seed = 1, k){
  # housekeeping
  function_out <- list()
  nchains <- length(fit)
  niter <- nrow(fit[[1]])
  function_out$niter <- niter
  function_out$nchains <- nchains
  
  # label-switching algorithm first
  samples_rl <- ord_ls(mcmc = fit, d = 2, seed = 3, force_K = k)
  
  # pseudo-convergence diagnostics
  out <- list()
  for(chain in 1:nchains){
    out[[chain]] <- with(samples_rl, {
      which_iter[
        (which_iter <= niter*chain) &
          (which_iter >= niter*(chain-1))
      ]
    })
  }
  ndx <- sapply(out, function(x) x[1:min(sapply(out, length))])
  if(any(sapply(out, length) < 100)) return(NA)
  ndx <- ndx[-c(1, nrow(ndx)),]
  all_samples <- cbind(samples_rl$mu, samples_rl$other)
  all_samples <- all_samples[,-which(colnames(all_samples) == "theta[1, 2]")]
  conv_diag <- matrix(NA, nrow = ncol(all_samples), ncol = 3)
  for(col in 1:ncol(all_samples)){
    tmp <- all_samples[which(samples_rl$which_iter %in% c(ndx)), col]
    tmp2 <- matrix(tmp, nrow = nrow(ndx))
    conv_diag[col, 1:3] <- c(
      Rhat(tmp2), ess_bulk(tmp2), ess_tail(tmp2)
    )
  }
  conv_diag_tbl <- tibble(
    param = colnames(all_samples),
    rhat = conv_diag[,1],
    ess_bulk = conv_diag[,2],
    ess_tail = conv_diag[,3]
  )
  function_out$conv <- conv_diag_tbl
  
  # number of clusters based on raw mcmc
  clus_id <- lapply(fit, function(x) x[,which(grepl("clus_id", colnames(x)))])
  clus_id_all <- do.call("rbind", clus_id)
  ngroups_mcmc <- apply(clus_id_all, 1, function(x) length(unique(x)))
  tab <- round(table(ngroups_mcmc) / length(ngroups_mcmc), 3)
  function_out$raw_nclus <- as.numeric(names(which(tab == max(tab))))
  
  # modal cluster assignments and means for z's
  z_mcmc <- samples_rl$other[,which(grepl("z", colnames(samples_rl$other)))]
  z_mcmc_dim1 <- z_mcmc[,grepl(", 1", colnames(z_mcmc))]
  z_mcmc_dim2 <- z_mcmc[,grepl(", 2", colnames(z_mcmc))]
  z_tbl <- tibble(
    obs_id = 1:ncol(z_mcmc_dim1),
    z1 = colMeans(z_mcmc_dim1),
    z2 = colMeans(z_mcmc_dim2),
    modal_clus = apply(samples_rl$clus_id, 2, getmode),
    clus_prob = colMeans(
      samples_rl$clus_id == matrix(rep(apply(samples_rl$clus_id, 2, getmode), each = nrow(samples_rl$clus_id)), nrow = nrow(samples_rl$clus_id))
    )
  )
  function_out$z_tbl <- z_tbl
  function_out$modal_nclus <- length(unique(z_tbl$modal_clus))
  function_out$nrow_relabeled <- length(samples_rl$which_iter)
  
  # out
  return(function_out)
}
prop_correct_pairwise <- function(true_clusters, est_clusters){
  true <- as.matrix(dist(true_clusters)) # 15 from data gen
  true[true>0] <- -1
  modeled <- as.matrix(dist(est_clusters))
  modeled[modeled>0] <- -1
  mean(true[upper.tri(true)] == modeled[upper.tri(modeled)])
}
get_conv_k1 <- function(coral_fit_samples){
  conv_diag <- matrix(NA, nrow = ncol(coral_fit_samples[[1]]), ncol = 3)
  for(col in 1:ncol(coral_fit_samples[[1]])){
    tmp <- matrix(NA, nrow(coral_fit_samples[[1]]), length(coral_fit_samples))
    for(i in 1:length(coral_fit_samples)){
      tmp[,i] <- coral_fit_samples[[i]][,col] 
    }
    conv_diag[col, 1:3] <- c(
      Rhat(tmp), ess_bulk(tmp), ess_tail(tmp)
    )
  }
  conv_diag_tbl <- tibble(
    param = colnames(coral_fit_samples[[1]]),
    rhat = conv_diag[,1],
    ess_bulk = conv_diag[,2],
    ess_tail = conv_diag[,3]
  )
  return(conv_diag_tbl)
}
get_ztbl_k1 <- function(coral_fit_samples){
  all_samples <- do.call("rbind", coral_fit_samples)
  z_samples <- all_samples[,grepl("z[[]", colnames(all_samples))]
  dim1 <- z_samples[,grepl(", 1", colnames(z_samples))]
  dim2 <- z_samples[,grepl(", 2", colnames(z_samples))]
  out <- tibble(
    obs_id = 1:ncol(dim1),
    z1 = colMeans(dim1),
    z2 = colMeans(dim2)
  ) %>%
    mutate(modal_clus = 1, clus_prob = 1)
  return(out)
}
get_waic_k1 <- function(Y, coral_fit, verbose = F){
  # combine chains
  all_samples <- do.call("rbind", coral_fit)
  
  # get params
  phi <- all_samples[,grepl("phi", colnames(all_samples))]
  alpha <- all_samples[,grepl("alpha", colnames(all_samples))]
  beta <- all_samples[,grepl("beta", colnames(all_samples))]
  z <- all_samples[,grepl("z", colnames(all_samples))]
  theta <- all_samples[,grepl("theta", colnames(all_samples))]
  
  # get dimension
  d <- sapply(colnames(z), function(x){
    str_sub(x, start = nchar(x)-1, end = nchar(x)-1)
  }) %>% as.numeric %>% max
  
  # compute lppd
  ## lik = alpha[i] + beta[j] + z[i] theta[j] + gamma[i]
  log_py <- matrix(NA, nrow = nrow(z), ncol = length(c(Y)))
  if(verbose) message("Calculating log likelihood")
  for(iter in 1:nrow(log_py)){
    alpha_vec_s <- alpha[iter,]
    beta_vec_s <- beta[iter,]
    z_matrix_s <- matrix(z[iter,], ncol = d)
    theta_matrix_s <- matrix(theta[iter,], ncol = d)
    phi_s <- phi[iter]
    
    linpred <- rep(alpha_vec_s, ncol(Y)) +
      rep(beta_vec_s, each = nrow(Y)) +
      c(z_matrix_s %*% t(theta_matrix_s)) 
    lambda <- unname(exp(linpred))
    nbp <- phi_s / (phi_s + lambda)
    
    # fix numerical issues
    phi_s[which(phi_s > 0.999)] <- 0.999
    phi_s[which(phi_s < 0.001)] <- 0.001
    
    log_py[iter,] <- phi_s * log(nbp) + c(Y) * log(1 - nbp)
  }
  
  # compute lppd
  if(verbose) message("Computing lppd")
  
  lppd <- sum(log(colMeans(exp(log_py))))
  
  # compute penalty
  if(verbose) message("Computing penalty")
  pwaic2 <- sum(apply(log_py, 2, function(x) var(x)))
  
  # compute waic
  waic <- -2*(lppd - pwaic2)
  out <- waic
  return(out)
}

# functions for simulation
# mat = compas_data[[1]][[1]]$Y
# true_coords = compas_data[[1]][[1]]$true_coords
# k = 2
# sim = 1
# niter = 100000
# nburnin = 50000
# thin = 5
# nchains = 3
# max_clus = 8
# max_attempts = 6
one_sim <- function(mat, true_coords, k, sim, niter = 75000, nchains = 3, nburnin = 25000, thin = 5, max_clus = 8, max_attempts = 10){
  # function to run one iteration of simulation
  out <- list()
  procrust <- matrix(NA, 3, 2)
  pairwise_clus_prop <- matrix(NA, 6, 1)
  nclus_out <- matrix(NA, 6, 1)
  conv_check <- matrix(1, 3, 1)
  
  # start with ordination techniques: dpord, coral, nmds
  ## dpord 
  ### "automate" model conv - only do this if looping through
  ### a bunch of simulated data sets
  dpord_rhat <- 999
  nattempts <- 1
  message("Fitting dpord")
  while(dpord_rhat > 1.1 & nattempts < max_attempts){
    message(paste0("Attempt "), nattempts)
    dpord_start <- Sys.time()
    this_cluster <- makeCluster(nchains)
    dpord_fit <- parLapply(
      cl = this_cluster,
      X = 1:nchains + 100*nattempts,
      fun = fit_model,
      code = dpord_code,
      data = list(
        Y = mat
      ),
      constants = list(
        nsites = nrow(mat),
        nspecies = ncol(mat),
        d = 2,
        max_clus = nrow(mat),
        S = diag(2),
        mu0 = rep(0, 2),
        Lambda0 = diag(2)
      ),
      inits = init_func(mat, d = 2, max_clus),
      niter = niter,
      burnin = nburnin,
      nchains = 1,
      thin = thin
    )
    stopCluster(this_cluster)
    dpord_end <- Sys.time()
    dpord_sum <- try(
      summarize_dpord(dpord_fit), 
      silent = T
    )
    if(length(dpord_sum) == 1){
      dpord_rhat <- 999
    } else{
      dpord_sum$runtime <- as.numeric(dpord_end - dpord_start, units = "secs")
      dpord_rhat <- max(dpord_sum$conv$rhat, na.rm = T)
    }
    saveRDS(dpord_fit, paste0( "simulations/models_raw/dpord_fit_", k, "_", sim, ".rds"))
    saveRDS(dpord_sum, paste0("simulations/summaries/dpord_sum_", k, "_", sim, ".rds"))
    nattempts <- nattempts + 1
  }

  # summaries
  if(length(dpord_sum) == 1){
    procrust[1,1] <- NA
    procrust[1,2] <- NA
    pairwise_clus_prop[1,1] <- NA
    nclus_out[1,1] <- NA
    conv_check[1,] <- 0
  } else{
    tmp <- vegan::procrustes(
      X = true_coords, Y = with(dpord_sum$z_tbl, cbind(z1, z2)),
      symmetric = TRUE
    )
    procrust[1,1] <- summary(tmp)$ss
    tmp <- vegan::procrustes(
      X = mat, Y = with(dpord_sum$z_tbl, cbind(z1, z2)),
      symmetric = TRUE
    )
    procrust[1,2] <- summary(tmp)$ss
    pairwise_clus_prop[1,1] <- prop_correct_pairwise(
      rep(1:k, each = 15),
      dpord_sum$z_tbl$modal_clus
    )
    nclus_out[1,1] <- dpord_sum$modal_nclus
  }
  
  ## coral model (fit k = 1, ..., max_clus)
  ### using nimble code I wrote for fair assessment
  waic <- matrix(NA, max_clus, 1)
  runtimes <- matrix(NA, max_clus, 1)
  for(ndx in 1:max_clus){
    message("Fitting coral with K = ", ndx)
    
    if(ndx == 1){
      coral_rhat <- 2
      nattempts <- 1
      while(coral_rhat > 1.1 & nattempts < max_attempts){
        message(paste0("Attempt "), nattempts)
        coral_start <- Sys.time()
        this_cluster <- makeCluster(nchains)
        coral_fit <- parLapply(
          cl = this_cluster,
          X = 1:nchains + nattempts*100,
          fun = fit_model_k1,
          code = coral_code_k1,
          data = list(
            Y = mat
          ),
          constants = list(
            nsites = nrow(mat),
            nspecies = ncol(mat),
            d = 2,
            S = diag(2),
            mu0 = rep(0, 2)
          ),
          inits = init_func(mat, d = 2, max_clus = 1),
          niter = niter,
          burnin = nburnin,
          nchains = 1,
          thin = thin
        )
        stopCluster(this_cluster)
        coral_end <- Sys.time()
        
        # separate waic and samples
        coral_sum <- list(
          niter = nrow(coral_fit[[1]]),
          nchains = length(coral_fit),
          conv = get_conv_k1(coral_fit),
          nclus = 1,
          runtime = as.numeric(coral_end - coral_start, units = "secs"),
          z_tbl = get_ztbl_k1(coral_fit),
          waic = get_waic_k1(mat, coral_fit, verbose = F)
        )
        
        coral_rhat <- max(coral_sum$conv$rhat, na.rm = T)
        nattempts <- nattempts + 1
        saveRDS(coral_fit, paste0("simulations/models_raw/coral_fit", ndx, "_", k, "_", sim, ".rds"))
        saveRDS(coral_sum, paste0("simulations/summaries/coral_sum", ndx, "_", k, "_", sim, ".rds"))
      }
      waic[ndx,1] <- coral_sum$waic
      runtimes[ndx,1] <- coral_sum$runtime
    } else{
      coral_rhat <- 2
      nattempts <- 1
      while(coral_rhat > 1.1 & nattempts < max_attempts){
        message(paste0("Attempt "), nattempts)
        coral_start <- Sys.time()
        this_cluster <- makeCluster(nchains)
        coral_fit <- parLapply(
          cl = this_cluster,
          X = 1:nchains + nattempts*100,
          fun = fit_model,
          code = coral_code,
          data = list(
            Y = mat
          ),
          constants = list(
            nsites = nrow(mat),
            nspecies = ncol(mat),
            d = 2,
            max_clus = nrow(mat),
            K = ndx,
            S = diag(2),
            mu0 = rep(0, 2),
            Lambda0 = diag(2),
            ones = rep(1, ndx)
          ),
          inits = init_func(mat, d = 2, ndx),
          niter = niter,
          burnin = nburnin,
          nchains = 1,
          thin = thin
        )
        stopCluster(this_cluster)
        coral_end <- Sys.time()
        
        coral_sum <- try(
          summarize_coral(coral_fit, k = ndx), 
          silent = T
        )
        if(length(coral_sum) == 1){
          waic[ndx,1] <- NA
          runtimes[ndx,1] <- as.numeric(coral_end - coral_start, units = "secs")
        } else{
          coral_sum$runtime <- as.numeric(coral_end - coral_start, units = "secs")
          coral_sum$waic <- get_waic_k1(mat, coral_fit)
          coral_rhat <- max(coral_sum$conv$rhat, na.rm = T)
          waic[ndx,1] <- coral_sum$waic
          runtimes[ndx,1] <- coral_sum$runtime
        }
        nattempts <- nattempts + 1
        saveRDS(coral_fit, paste0("simulations/models_raw/coral_fit", ndx, "_", k, "_", sim, ".rds"))
        saveRDS(coral_sum, paste0("simulations/summaries/coral_sum", ndx, "_", k, "_", sim, ".rds"))
      }
    }
  }
  
  # select "best" model
  if(sum(!is.na(waic)) == 0){
    procrust[2,1] <- NA
    procrust[2,2] <- NA
    pairwise_clus_prop[2,1] <- NA
    nclus_out[2,1] <- NA
    conv_check[2,1] <- 0
  } else if(sum(!is.na(waic[-1,])) == 0){
    conv_check[2,1] <- 0
    nclus_out[2,1] <- 1
    waic_ <- na.omit(waic)
    coral_sum <- readRDS(
      paste0(
        "simulations/summaries/coral_sum", 
        1,
        "_", k, "_", sim, ".rds"
      )
    )
    
    # summaries
    ## procrustes
    tmp <- vegan::procrustes(
      X = true_coords, Y = with(coral_sum$z_tbl, cbind(z1, z2)),
      symmetric = TRUE
    )
    procrust[2,1] <- summary(tmp)$ss
    
    tmp <- vegan::procrustes(
      X = mat, Y = with(coral_sum$z_tbl, cbind(z1, z2)),
      symmetric = TRUE
    )
    procrust[2,2] <- summary(tmp)$ss
    
    ## prop of correct pairwise group_membership
    pairwise_clus_prop[2,1] <- prop_correct_pairwise(
      rep(1:k, each = 15),
      coral_sum$z_tbl$modal_clus
    )
  } else{
    coral_sum <- readRDS(
      paste0(
        "simulations/summaries/coral_sum", 
        which(c(waic) == min(c(waic), na.rm = T)),
        "_", k, "_", sim, ".rds"
      )
    )
    
    # summaries
    ## procrustes
    tmp <- vegan::procrustes(
      X = true_coords, Y = with(coral_sum$z_tbl, cbind(z1, z2)),
      symmetric = TRUE
    )
    procrust[2,1] <- summary(tmp)$ss
    
    tmp <- vegan::procrustes(
      X = mat, Y = with(coral_sum$z_tbl, cbind(z1, z2)),
      symmetric = TRUE
    )
    procrust[2,2] <- summary(tmp)$ss
    
    ## prop of correct pairwise group_membership
    pairwise_clus_prop[2,1] <- prop_correct_pairwise(
      rep(1:k, each = 15),
      coral_sum$z_tbl$modal_clus
    )
    nclus_out[2,1] <- coral_sum$modal_nclus
  }
  
  # nmds
  message("Running nmds")
  nmds_start <- Sys.time()
  capture.output(nmds <- vegan::metaMDS(mat), file = "NUL")
  nmds_end <- Sys.time()
  nmds_sum <- list(
    z_tbl = tibble(
      obs_id = nrow(nmds$points),
      z1 = nmds$points[,1],
      z2 = nmds$points[,2]
    ),
    runtime = as.numeric(nmds_end - nmds_start, units = "secs")
  )
  conv_check[3,] <- ifelse(
    nmds$maxits == nmds$iters, 0, 1
  )
  saveRDS(nmds, paste0( "simulations/models_raw/nmds_", k, "_", sim, ".rds"))
  saveRDS(nmds_sum, paste0("simulations/summaries/nmds_sum_", k, "_", sim, ".rds"))
  
  tmp <- vegan::procrustes(
    X = true_coords, Y = nmds$points,
    symmetric = TRUE
  )
  procrust[3,1] <- summary(tmp)$ss
  
  tmp <- vegan::procrustes(
    X = mat, Y = nmds$points,
    symmetric = TRUE
  )
  procrust[3,2] <- summary(tmp)$ss
  
  # summarize ordination techniques
  ordination <- tibble(
    method = c('dpord', "coral", "nmds"),
    procrustes_true = c(procrust[,1]),
    procrustes_data = c(procrust[,2]),
    runtime = c(as.numeric(dpord_end - dpord_start, units = "secs"), sum(c(runtimes)), nmds_sum$runtime),
    conv = c(conv_check)
  ) %>%
    mutate(
      k = k,
      sim = sim
    )
  saveRDS(
    ordination, 
    file = paste0(
      "simulations/agg_summaries/ordination_", k, "_", "sim", ".rds"
    )
  )
  
  # clustering techniques
  # algorithmic techniques
  ## maximize silhouette distance with ward clustering
  bray_dist <- as.matrix(vegan::vegdist(mat, method = "bray"))
  ward <- cluster::agnes(
    x = bray_dist,
    diss = TRUE,
    method = "ward"
  )
  ward_sil <- sapply(
    2:max_clus, function(x){
      mean(silhouette(cutree(ward, k = x), bray_dist)[,3])
    }
  )
  ward_nclus <- (2:max_clus)[which(ward_sil == max(ward_sil))]
  ward_clus <- cutree(ward, ward_nclus)
  pairwise_clus_prop[3,1] <- prop_correct_pairwise(
    rep(1:k, each = 15),
    ward_clus
  )
  nclus_out[3,1] <- ward_nclus

  ## pam with bray curtis, maximizing sil dist
  pam <- lapply(
    2:max_clus, function(x){
      cluster::pam(
        x = bray_dist,
        k = x,
        diss = "TRUE",
      )
    }
  )
  pam_sil <- sapply(
    pam, function(x){
      mean(silhouette(x$clustering, bray_dist)[,3])
    }
  )
  pam_nclus <- (2:max_clus)[which(pam_sil == max(pam_sil))]
  pam_clus <- pam[[which(pam_sil == max(pam_sil))]]$clustering
  pairwise_clus_prop[4,1] <- prop_correct_pairwise(
    rep(1:k, each = 15),
    pam_clus
  )
  nclus_out[4,1] <- pam_nclus
  
  ## k-means with cal
  ch <- sapply(
    2:max_clus, function(x){
      kmeans <- kmeans(mat, x)
      vegan::cIndexKM(kmeans, mat)[1]
    }
  )
  ssi <- sapply(
    2:max_clus, function(x){
      kmeans <- kmeans(mat, x)
      vegan::cIndexKM(kmeans, mat)[2]
    }
  )
  
  pairwise_clus_prop[5, 1] <- prop_correct_pairwise(
    rep(1:k, each = 15),
    kmeans(mat, (2:max_clus)[which(ch == max(ch))])$cluster
  )
  nclus_out[5,1] <- (2:max_clus)[which(ch == max(ch))]
  
  pairwise_clus_prop[6, 1] <- prop_correct_pairwise(
    rep(1:k, each = 15),
    kmeans(mat, (2:max_clus)[which(ssi == max(ssi))])$cluster
  )
  nclus_out[6,1] <- (2:max_clus)[which(ssi == max(ssi))]
  
  ordination <- tibble(
    method = c('dpord', "coral", "ward", "pam", "ch", "ssi"),
    pairwise_clus_prop = c(pairwise_clus_prop),
    nclus = c(nclus_out)
  ) %>%
    mutate(
      k = k,
      sim = sim
    )
  saveRDS(
    clustering, 
    file = paste0(
      "simulations/agg_summaries/clustering_", k, "_", "sim", ".rds"
    )
  )
  
  return(
    list(
      ordination = ordination,
      clustering = clustering
    )
  )
  
}
