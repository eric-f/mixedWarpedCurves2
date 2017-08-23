rm(list=ls())
gc()

library(mixedWarpedCurves2)
library(plyr)
library(ggplot2)

# dat0 <- readRDS("../data/test_data.rds")
dat0 <- readRDS("data/test_data.rds")
# dat0$y <- dat0$y * attr(dat0$y, "scaled:scale") + attr(dat0$y, "scaled:center")

library(SESfda)
dat0 <- genDataWiggle(n = 20,
                      ni = 1000,
                      basis.control = control.bspline.basis(f.nknots=2),
                      Sigma0 = diag(c(20^2, 0.0025)))
dat0$y <- dat0$y.shape + dat0$y.noise
ggplot(dat0) +
  geom_line(aes(x=x, y=y, col=as.factor(id)))

# dat0 <- readRDS("~/org/project/Simulations/data/w1-20x100/data_w1-20x100-set_000.rds")
my_fit <-
  mixedWarpedCurves2::fsim(y = dat0$y,
                           obs_time = dat0$x,
                           curve_id = dat0$id,
                           n_clust = 1,
                           saem_control = control_saem(n_saem_iter = 500,
                                                       n_saem_burn = 100,
                                                       saem_step_seq_pow = 0.75,
                                                       n_mcmc_burn = 5,
                                                       n_core = 1,
                                                       accept_rate_window = 5,
                                                       # prop_sigma = 7.8125e-05,
                                                       prop_sigma = 1e-2,
                                                       need_centering = TRUE,
                                                       accept_rate_lb = 0.17,
                                                       accept_rate_ub = 0.33,
                                                       h_n_knots = 3+2,
                                                       f_n_knots = 2+2,
                                                       ind_amp = TRUE))

round(
  cbind(
  "Truth" = c(attr(dat0, "par0")$alp, attr(dat0, "par0")$sigma2, diag(attr(dat0, "par0")$Sigma), attr(dat0, "par0")$tau),
  "Estimate" = c(my_fit$pars$alpha, my_fit$pars$sigma2, diag(my_fit$pars$big_sigma), my_fit$pars$tau_clusters),
  "SE" = sqrt(diag(my_fit$se_info$varcov))), 3)

diag(my_fit$se_info$scaling_mat)

my_fitted <- ldply(my_fit$curve, function(crv){
  data.frame("x"=crv$x,
             "y"=crv$y,
             "warped_x"=crv$warped_x,
             "fitted_y"=crv$fitted_y,
             "curve_id"=crv$curve_id)
})

ggplot(my_fitted) +
  geom_line(aes(x=x, y=y, col=as.factor(curve_id))) +
  geom_line(aes(x=x, y=fitted_y, group=curve_id), linetype=2) +
  facet_wrap(~curve_id)


sigma_a_tmp <- sapply(my_fit$curves, function(x){x$sapprox_a})
plot(t(sigma_a_tmp), xlim=c(-60, 60), ylim=c(0.8, 1.2))
points(attributes(dat0)$RE[,1:2], col="red")

plot(sigma_a_tmp[1,], attributes(dat0)$RE[,1])
plot(sigma_a_tmp[2,], attributes(dat0)$RE[,2])
