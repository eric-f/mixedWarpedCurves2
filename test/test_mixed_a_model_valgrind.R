rm(list=ls())
gc()

# source("../R/sim_warping_mixture.R")
library(mixedWarpedCurves2)
library(ggplot2)
library(gtools)

kappa0 <- c(1,2,2,1)
kappa0 <- kappa0 / sum(kappa0)
kappa1 <- c(0,1,1,0)
kappa1 <- kappa1 / sum(kappa1)
kappa2 <- c(1.2,1,1,1.2)
kappa2 <- kappa2 / sum(kappa2)
kappa3 <- c(0,1,2,1)
kappa3 <- kappa3 / sum(kappa3)
kappa4 <- c(1,1,0,0)
kappa4 <- kappa4 / sum(kappa4)

dat0 <- sim_warping_mixture(200, rep(1/3, 3),
                            rbind(kappa1,
                                  kappa2,
                                  kappa3),
                            ni = 1000,
                            tau = 5,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)

library(lineprof)

saem_control <- control_saem(n_saem_iter = 50,
                             n_saem_burn = 20,
                             saem_step_seq_pow = 1,
                             n_mcmc_burn = 1,
                             n_core = 1,
                             accept_rate_window = 5,
                             prop_sigma = 1e-2,
                             need_centering = FALSE,
                             accept_rate_lb = 0.17,
                             accept_rate_ub = 0.33,
                             f_n_knots = 1+5)

## Random start
# dat0$clust <- rep(sample(0:2, 200, replace=TRUE), each=201)
# tmp <- lineprof::lineprof({
# init_clust <- rep(sample(3, 200, replace=TRUE), 1000) - 1
init_clust = as.integer(dat0$clust) - 1
system.time(my_fit <-
              mixedWarpedCurves2::fsim_mixture_a_model(
                y = dat0$y,
                obs_time = dat0$x,
                curve_id = dat0$id,
                init_clust = init_clust,
                n_clust = 3,
                saem_control = saem_control,
                trace = FALSE))
# })
# tmp
# lineprof::shine(tmp)

matplot(t(my_fit$pars_track$sampled_m_track), type="l")

pred_cluster <- sapply(my_fit$curve, function(crv){
  which.max(crv$sapprox_cluster_membership)
  })
true_cluster <- as.integer(dat0$clust[dat0$x==0])

table(pred_cluster, true_cluster)

max(apply(apply(apply(my_fit$pars_track$sampled_m_track, 1, diff), 2, abs), 2, which.max))
