rm(list=ls())
gc()
library(mixedWarpedCurves2)
library(plyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(gtools)


kappa1 <- c(0,1,1,0)
kappa1 <- kappa1 / sum(kappa1)
kappa2 <- c(1.2,1,1,1.2)
kappa2 <- kappa2 / sum(kappa2)
kappa3 <- c(0,1,2,1)
kappa3 <- kappa3 / sum(kappa3)
dat0 <- sim_warping_mixture(200, rep(1/3, 3),
                            rbind(kappa1,
                                  kappa2,
                                  kappa3),
                            ni = 101,
                            tau = 10,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)
init_clust <- as.integer(dat0$clust[dat0$x==0]) - 1
y <- matrix(dat0$y, ncol=200)
x <- dat0$x[1:101]

kmax <- 10

llks <- numeric(kmax)
is_out_lst <- mod_fit_lst <- rep(list(NULL), kmax)

for(k in seq(along=llks)){
  print(k)
  if(k!=1){
    bs_out <- mixedWarpedCurves2::mixed_shape(y, x, nCat = k, nTry = 30)
    bs_clust <- apply(bs_out$p_jk, 1, which.max) - 1
  }
  else{
    bs_clust <- rep(0, 200)
  }
  system.time(my_fit <-
                mixedWarpedCurves2::fsim_unimodal(
                  y = dat0$y,
                  obs_time = dat0$x,
                  curve_id = dat0$id,
                  init_clust = bs_clust,
                  n_clust = k,
                  saem_control = control_saem(n_saem_iter = 1000,
                                              n_saem_burn = 100,
                                              saem_step_seq_pow = 1,
                                              n_mcmc_burn = 5,
                                              n_core = 1,
                                              accept_rate_window = 5,
                                              prop_sigma = 1e-2,
                                              need_centering = FALSE,
                                              accept_rate_lb = 0.17,
                                              accept_rate_ub = 0.33,
                                              h_n_knots = 1+4),
                  trace = FALSE))
  mod_fit_lst[[k]] <- my_fit
}

seed <- round(runif(1)*1e6)
for(k in seq(along=llks)){
  print(k)
  is_out_lst[[k]] <- is_out <- logLik_fsim(my_fit, n_sim=100000,
                                           method = "joint_mc", seed = seed)
  llks[k] <- sum(sapply(is_out$curves, function(x){
      log(mean(exp(x$mc_log_weight + x$mc_logLik)))}))
  # llks[k] <- sum(sapply(is_out$curves, function(x){
  #   scl <- max(is_out$curves[[1]]$mc_logLik)
  #   log(sum(exp(is_out$curves[[1]]$mc_logLik - scl - log(10000)))) + scl
  # }))
}
plot(llks)

w <- exp(is_out_lst[[1]]$curves[[6]]$mc_log_weight)
10000 * mean(w)^2 / mean(w^2)

