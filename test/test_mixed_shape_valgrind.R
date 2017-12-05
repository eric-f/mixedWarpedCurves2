# rm(list=ls())
# gc()

library(mixedWarpedCurves2)

dat0 <- readRDS("../data/test_data_mixed_warped_unimodal.rds")

## Random start
# init_clust <- rep(sample(3, 200, replace=TRUE), 1000) - 1
init_clust = as.integer(dat0$clust) - 1
system.time({
  my_fit <-
    mixedWarpedCurves2::fsim_mixed_shape(
      y = dat0$y,
      obs_time = dat0$x,
      curve_id = dat0$id,
      init_clust = init_clust,
      n_clust = 3,
      saem_control = control_saem(n_saem_iter = 1000,
                                  n_saem_burn = 100,
                                  saem_step_seq_pow = 1,
                                  n_mcmc_burn = 1,
                                  n_core = 1,
                                  accept_rate_window = 5,
                                  prop_sigma = 1e-2,
                                  need_centering = FALSE,
                                  accept_rate_lb = 0.17,
                                  accept_rate_ub = 0.33,
                                  f_n_knots = 1+5),
      trace = FALSE)
})
