rm(list=ls())
gc()
library(mixedWarpedCurves2)
library(plyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(gtools)

out <- readRDS("~/org/project/Simulations/module/sim-20171113_unimodal_mixture_brillinger_start_3-200x1000/output/unimodal_mixture_brillinger_start_3-200x1000_set_000.rds")

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
                            ni = 11,
                            tau = 10,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)
init_clust <- as.integer(dat0$clust[dat0$x==0]) - 1
system.time(my_fit <-
              mixedWarpedCurves2::fsim_unimodal(
                y = dat0$y,
                obs_time = dat0$x,
                curve_id = dat0$id,
                init_clust = init_clust,
                n_clust = 3,
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

saemObj <- my_fit
saemObj <- out$saemOut
n_mc = 1000

is_out <- logLik_fsim(saemObj, n_sim=n_mc, method = "joint_is", seed = 100)
range(exp(is_out$curves[[1]]$mc_log_weight))
range(is_out$curves[[1]]$mc_logLik)
plot(is_out$curves[[1]]$mc_rss)
log(mean(exp(is_out$curves[[1]]$mc_log_weight + is_out$curves[[1]]$mc_logLik)))
# range(is_out$curves[[1]]$mc_logLik)
# is_out$curves[[1]]$mc_log_weight
# (is_out$curves[[1]]$mc_log_weight_a + log(colSums(exp(is_out$curves[[1]]$mc_log_weight_w))))

range(is_out$curves[[1]]$mc_log_weight_a)
range(is_out$curves[[1]]$mc_log_weight_w)
plot(c(is_out$curves[[1]]$mc_log_weight),
       c(is_out$curves[[1]]$mc_log_weight_w[1,]))
plot(c(is_out$curves[[1]]$mc_log_weight),
     c(is_out$curves[[1]]$mc_log_weight_w[2,]))
plot(c(is_out$curves[[1]]$mc_log_weight),
     c(is_out$curves[[1]]$mc_log_weight_w[3,]))

w = exp(is_out$curves[[1]]$mc_log_weight)
n_mc * mean(w)^2 / mean(w^2)

hist(is_out$curves[[1]]$mc_log_weight)

log(mean(exp(is_out$curves[[1]]$mc_log_weight + is_out$curves[[1]]$mc_logLik)))



mc_out <- logLik_fsim(saemObj, n_sim=n_mc, method = "joint_mc")
plot(mc_out$curves[[1]]$mc_rss)
scl <- max(mc_out$curves[[1]]$mc_logLik)
log(sum(exp(mc_out$curves[[1]]$mc_logLik - scl - log(n_mc)))) + scl
# out$pars
# saemObj$pars
#
# str(out$curves[[1]])
# str(saemObj$curves[[1]])

# names(llk_out)
# llk_out$aux
# str(llk_out$curves[[1]])
#
# str(llk_out$curves[[1]])
# range(exp(llk_out$curves[[1]]$mc_sample))
#
# llk_out$curves[[2]]$w
# llk_out$pars[["kappa_clusters"]]

mcmc_out <- logLik_fsim(saemObj, n_sim=n_mc, method = "mcmc")
plot(mcmc_out$curves[[1]]$mc_rss)
range(mcmc_out$curves[[1]]$mc_logLik)
plot(mcmc_out$curves[[1]]$mc_logLik)
log(mean(exp(mcmc_out$curves[[1]]$mc_logLik[-(1:100)])))

-sum((saemObj$curves[[1]]$y - saemObj$curves[[1]]$fitted_y)^2)/saemObj$pars$sigma2/2

a=0.0001;curve(dnorm(x, sd=a, log = TRUE), from=-a, a)
plot(mcmc_out$curves[[1]]$a[1,], type="l");abline(h=mcmc_out$curves[[1]]$sapprox_a[1], col="red")
plot(mcmc_out$curves[[1]]$a[2,], type="l");abline(h=mcmc_out$curves[[1]]$sapprox_a[2], col="red")
plot(mcmc_out$curves[[1]]$w[2,], type="l");abline(h=mcmc_out$curves[[1]]$sapprox_w[2], col="red")
plot(mcmc_out$curves[[1]]$y~mcmc_out$curves[[1]]$x)
