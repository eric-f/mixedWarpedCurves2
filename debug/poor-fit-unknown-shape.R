library(mixedWarpedCurves2)
library(gtools)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

data_generator_beta_shape_mixture <- function(seed, k, n, ni, shape_a, shape_b){
  ## 1. Simulate wiggly data
  kappa0 <- c(1,2,2,1)
  kappa0 <- kappa0 / sum(kappa0)
  kappa1 <- c(0,1,1,0)
  kappa1 <- kappa1 / sum(kappa1)
  kappa2 <- c(1.2,1,1,1.2)
  kappa2 <- kappa2 / sum(kappa2)
  kappa3 <- c(0,1,2,1)
  kappa3 <- kappa3 / sum(kappa3)
  kappa4 <- c(1,2,1,0)
  kappa4 <- kappa4 / sum(kappa4)
  kappas = rbind(kappa1, kappa2, kappa3, kappa4)

  set.seed(seed)
  dat0 <- sim_warping_mixture(n, rep(1/k, k),
                              kappas[1:k,],
                              ni = ni,
                              tau = 20,
                              mu_sh = -25, mu_sc = 500,
                              sd_sh = 50, sd_sc=50, sd_err = 10,
                              shape_a = shape_a, shape_b = shape_b)
  dat0$seed <- seed
  return(dat0)
}

k <- 3
dat0 <- data_generator_beta_shape_mixture(seed=0, k=k, n=200, ni=100, 2, 2)
dat0$y <- dat0$y/ max(abs(dat0$y))

(h_knots0 <- seq(0,1,length=5))

x <- unique(dat0$x)
y <- dat0 %>%
  select(x, y, id) %>%
  spread(x, y) %>%
  remove_rownames %>%
  column_to_rownames("id") %>%
  t

seed <- 0
set.seed(seed)
bs97_fit <- mixed_shape(y=y, x=x, nCat=k, nTry=10)
init_clust <- apply(bs97_fit$p_jk, 1, which.max)

ggplot(dat0) +
  geom_line(aes(x=x, y=y, col=clust, group=id))

out <-
  mixedWarpedCurves2::fsim_mixed_warped_curves(
  # mixedWarpedCurves2::fsim_unimodal(
  # mixedWarpedCurves2::fsim_mixed_warping(
  # mixedWarpedCurves2::fsim(
    y = dat0$y,
    obs_time = dat0$x,
    curve_id = dat0$id,
    init_clust = init_clust,
    # init_clust = rep(1, 200),
    n_clust = k,
    saem_control = control_saem(
      seed = seed,
      n_saem_iter = 1000,
      n_saem_burn = 100,
      saem_step_seq_pow = 0.75,
      n_mcmc_burn = 5,
      n_core = 1,
      accept_rate_window = 5,
      prop_sigma = 1e-2,
      need_centering = FALSE,
      accept_rate_lb = 0.17,
      accept_rate_ub = 0.33,
      h_knots = h_knots0,
      ind_amp = TRUE))


fitt_dat <- plyr::ldply(out$curves, function(crv){
  data.frame(
    id = crv$curve_id,
    x=crv$x,
    y=crv$y,
    fitted_y=crv$fitted_y,
    clust=as.character(which.max(crv$sapprox_cluster_membership)))
})

fitt_dat %>%
  filter(id %in% sample(200,9)) %>%
  ggplot() +
  geom_line(aes(x=x, y=fitted_y, group=id), col="gold", size=2) +
  geom_line(aes(x=x, y=y, group=id, col=clust)) +
  facet_wrap(~id)

out$pars$sigma2

clust <- plyr::laply(out$curves, function(crv){which.max(crv$sapprox_cluster_membership)})

# clust_unimodal_h5 <- clust
# clust_unimodal_h3 <- clust
# clust_mxwp_h3 <- clust
# clust_mxwp_h5 <- clust
table(clust_unimodal_h3, clust_unimodal_h5)
table(clust_unimodal_h3, clust_mxwp_h3)
table(clust_unimodal_h5, clust_mxwp_h5)
table(clust_mxwp_h3, clust_mxwp_h5)

true_clust <- dat0$clust[dat0$x==0]

table(true_clust, init_clust)
table(true_clust, clust_unimodal_h3)
table(true_clust, clust_unimodal_h5)
table(true_clust, clust_mxwp_h3)
table(true_clust, clust_mxwp_h5)

mclust::adjustedRandIndex(true_clust, init_clust)
mclust::adjustedRandIndex(true_clust, clust_unimodal_h3)
mclust::adjustedRandIndex(true_clust, clust_unimodal_h5)
mclust::adjustedRandIndex(true_clust, clust_mxwp_h3)
mclust::adjustedRandIndex(true_clust, clust_mxwp_h5)

