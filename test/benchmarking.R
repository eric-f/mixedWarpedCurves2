rm(list=ls())
gc()

library(mixedWarpedCurves2)
library(mixedWarpedCurves)
library(microbenchmark)

# dat0 <- readRDS("../data/test_data.rds")
dat0 <- readRDS("data/test_data.rds")

microbenchmark(
  saemOut <- mixedWarpedCurves::fsim(
    y = dat0$y,
    obs_time = dat0$x,
    curve_id = dat0$id,
    init_pars = NULL,
    basis.control = attr(dat0, "basis.control"),
    pars.control = control_model_param(fixed.Sigma = FALSE),
    sa.control = control_sa(nIter = 100,
                            alphaSAEM = 0.75,
                            nCore = 1,
                            nBurnSAEM = 10,
                            nBurnMCMC = 5,
                            prop_sigma = 1e-2,
                            centering = TRUE,
                            accept_rate_lb = 0.17,
                            accept_rate_ub = 0.33),
    track_pars = FALSE),
  my_fit <-
    mixedWarpedCurves2::fsim(y = dat0$y,
                             obs_time = dat0$x,
                             curve_id = dat0$id,
                             n_clust = 1,
                             saem_control = control_saem(n_saem_iter = 100,
                                                         n_saem_burn = 5,
                                                         saem_step_seq_pow = 1,
                                                         n_mcmc_burn = 5,
                                                         n_core = 1,
                                                         accept_rate_window = 5,
                                                         # prop_sigma = 7.8125e-05,
                                                         prop_sigma = 1e-2,
                                                         need_centering = TRUE,
                                                         accept_rate_lb = 0.17,
                                                         accept_rate_ub = 0.33,
                                                         h_n_knots = 5+2,
                                                         f_n_knots = 5+2)),
  times=10L
)
