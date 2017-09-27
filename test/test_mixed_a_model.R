rm(list=ls())
gc()

source("R/sim_warping_mixture.R")
library(mixedWarpedCurves2)
library(ggplot2)

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
                            ni = 201,
                            tau = 20,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)

saem_control <- control_saem(n_saem_iter = 100,
             n_saem_burn = 50,
             saem_step_seq_pow = 1,
             n_mcmc_burn = 5,
             n_core = 1,
             accept_rate_window = 5,
             prop_sigma = 1e-2,
             need_centering = FALSE,
             accept_rate_lb = 0.17,
             accept_rate_ub = 0.33,
             f_n_knots = 1+5)

## Random start
# dat0$clust <- rep(sample(0:2, 200, replace=TRUE), each=201)
system.time(my_fit <-
              mixedWarpedCurves2::fsim_mixture_a_model(
                y = dat0$y,
                obs_time = dat0$x,
                curve_id = dat0$id,
                init_clust = as.integer(dat0$clust) - 1,
                n_clust = 3,
                saem_control = saem_control,
                trace = FALSE))

# my_fit$pars
# my_fit$curves[[1]]$fitted_y

help(perm.test)
my_fit$pars
f_knots <- sort(c(rep(range(saem_control$f_knots), saem_control$f_order-1), saem_control$f_knots))
x0 <- seq(0, 1, length=1001)
bfx <- splineDesign(f_knots, x0, saem_control$f_order, rep(0, length(x0)))
matplot(bfx %*% my_fit$pars$alpha, type="l")


fitted_mod <- ldply(my_fit$curves, function(crv){
  data.frame(fitted_y1=c(crv$post_amp_mu[1,1] + crv$post_amp_mu[2,1] *
                           crv$fitted_y[,1]),
             fitted_y2=c(crv$post_amp_mu[1,2] + crv$post_amp_mu[2,2] *
                           crv$fitted_y[,2]),
             fitted_y3=c(crv$post_amp_mu[1,3] + crv$post_amp_mu[2,3] *
                           crv$fitted_y[,3]),
             pred_clust=c(which.max(crv$sapprox_cluster_membership)))
  }) %>%
  mutate(
    fitted_y = fitted_y1 * (pred_clust==1) + fitted_y2 * (pred_clust==2) + fitted_y3 * (pred_clust==3),
    pred_clust = as.factor(pred_clust)
  )
dat1 <- cbind.data.frame(dat0, fitted_mod)

# ggplot(dat1) +
#   geom_line(aes(x=x, y=y, group=id, color=clust)) +
#   geom_line(aes(x=x, y=fitted_y, group=id, color=pred_clust), linetype=2) +
#   facet_wrap(~id)

ldply(my_fit$curves, function(crv){
  data.frame(id=crv$curve_id+1,
             pred_clust=c(which.max(crv$sapprox_cluster_membership)))
}) %>%
  join(., dat0 %>% dplyr::filter(x==0)) %>%
  xtabs(~clust + pred_clust, data=.)


matplot(t(my_fit$pars_track$alpha_track[,1,]), type="l")
matplot(t(my_fit$pars_track$alpha_track[,2,]), type="l")
matplot(t(my_fit$pars_track$alpha_track[,3,]), type="l")

x <- sample(6, 1000, replace = T)
y <- sample(6, 1000, replace = T)

perm.test(x, y, exact = TRUE)
