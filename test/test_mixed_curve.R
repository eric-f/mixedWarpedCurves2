# rm(list=ls())
# gc()

library(plyr)
library(dplyr)
library(gtools)
library(ggplot2)
library(mixedWarpedCurves2)

kappa0 <- c(1,2,2,1)
kappa0 <- kappa0 / sum(kappa0)
kappa1 <- c(0,1,1,0)
kappa1 <- kappa1 / sum(kappa1)
kappa2 <- c(1.5,1,1,1.5)
kappa2 <- kappa2 / sum(kappa2)
kappa3 <- c(0,1,2,1)
kappa3 <- kappa3 / sum(kappa3)
kappa4 <- c(1,1,0,0)
kappa4 <- kappa4 / sum(kappa4)

dat0 <- sim_warping_mixture(200, rep(1/3, 3),
                            rbind(kappa1,
                                  kappa2,
                                  kappa3),
                            ni = 201,
                            tau = 20,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)
dat0 %>%
  ggplot() +
  geom_line(aes(x=x, y=y, group=id, col=as.factor(clust)),
            show.legend = FALSE)


## Random start
init_clust = as.integer(dat0$clust) - 1
system.time({
  my_fit <-
    mixedWarpedCurves2::fsim_mixed_warping(
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
                                  f_n_knots = 9,
                                  h_n_knots = 3),
      trace = FALSE)
})



### Dev
y <- dat0$y
obs_time = dat0$x
curve_id = dat0$id
init_clust = as.integer(dat0$clust) - 1
n_clust = 3

saem_control = control_saem(n_saem_iter = 1000,
                            n_saem_burn = 100,
                            saem_step_seq_pow = 1,
                            n_mcmc_burn = 1,
                            n_core = 1,
                            f_n_knots = 9,
                            h_n_knots = 3)
