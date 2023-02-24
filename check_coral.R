# check coral
library(ade4)
data(doubs)
Y <- doubs$fish[-8,] %>% as.matrix
Y[which(Y>0)] <- 1

source("coral.R")

coral.fit <- fit.coral(
  y = Y,
  family = "binomial",
  num.lv = 2,
  num.clus = 2,
  site.eff = TRUE,
  save.model = TRUE,
  n.burnin = nburnin,
  n.iter = niter,
  n.thin = thin,
  seed = 1
)

plot(coral.fit$lv.mean, col = coral.fit$clus.modal, type = "n")
text(coral.fit$lv.mean, labels = 1:nrow(coral.fit$lv.mean), col = coral.fit$clus.modal, type = "n")
