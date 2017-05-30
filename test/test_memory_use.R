rm(list=ls())
gc()

library(mixedWarpedCurves2)
# library(mixedWarpedCurves)
library(ggplot2)
library(plyr)

dat0 <- readRDS("../data/test_data.rds")
# dat0 <- readRDS("data/test_data.rds")
# dat0 <- readRDS("~/org/lib/mixedWarpedCurves/data/test_data.rds")
# dat0 <- readRDS("~/org/project/Simulations/data/w1-20x1000/data_w1-20x1000-set_000.rds")
# system.time({
# # lineprof({
#   set.seed(1)
#   saemOut <- mixedWarpedCurves::fsim(
#     y = dat0$y,
#     obs_time = dat0$x,
#     curve_id = dat0$id,
#     init_pars = NULL,
#     basis.control = attr(dat0, "basis.control"),
#     pars.control = control_model_param(fixed.Sigma = FALSE),
#     sa.control = control_sa(nIter = 100,
#                             alphaSAEM = 0.75,
#                             nCore = 1,
#                             nBurnSAEM = 10,
#                             nBurnMCMC = 5,
#                             prop_sigma = 1e-2,
#                             centering = TRUE,
#                             accept_rate_lb = 0.17,
#                             accept_rate_ub = 0.33),
#     track_pars = FALSE)
# })
#
#
# saemOutDf <- predict_fsim(saemOut)
# saemOutDf$id <- factor(saemOutDf$id, levels = 1:20)
# ggplot(saemOutDf) +
#   geom_line(aes(x=x, y=y, col=id)) +
#   geom_line(aes(x=x, y=fitted_y, group=id), linetype=2) +
#   facet_wrap(~id)

bctrl <- attr(dat0, "basis.control")

system.time(my_fit <-
              mixedWarpedCurves2::fsim(y = dat0$y,
                                       obs_time = dat0$x,
                                       curve_id = dat0$id,
                                       saem_control = control_saem(n_saem_iter = 1000,
                                                                   n_saem_burn = 100,
                                                                   saem_step_seq_pow = 1,
                                                                   n_mcmc_burn = 5,
                                                                   n_core = 4,
                                                                   accept_rate_window = 5,
                                                                   # prop_sigma = 7.8125e-05,
                                                                   prop_sigma = 1e-2,
                                                                   need_centering = TRUE,
                                                                   accept_rate_lb = 0.17,
                                                                   accept_rate_ub = 0.33,
                                                                   h_n_knots = bctrl$h.nknots+2,
                                                                   f_n_knots = bctrl$f.nknots+2)))
# saveRDS(my_fit, file = "my_saem_fit.rds")
plot(my_fit$aux$proposal_sigma_history, type="l")
plot(my_fit$aux$mh_accept_rate_history)
plot(my_fit$pars_track$sigma2_track[-(1:80)], type="l")
plot(my_fit$pars_track$tau_track[-(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[1, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[2, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[3, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[4, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[5, -(1:80)], type="l")
matplot(t(my_fit$pars_track$big_sigma_track)[-(1:80),1], type="l")
matplot(t(my_fit$pars_track$big_sigma_track)[-(1:80),2], type="l")
matplot(t(my_fit$pars_track$big_sigma_track)[-(1:80),4], type="l")
out <- ldply(my_fit$curves, function(x){
  data.frame(id=x$curve_id,
             x=x$x,
             y=x$y,
             warped_x=x$warped_x,
             fitted_y=x$fitted_y)
})
out$id <- as.factor(out$id+1)
#
# var(dat0$y - dat0$y.shape)
# var(saemOutDf$y * sd(dat0$y) + mean(dat0$y) - saemOutDf$fitted_y * sd(dat0$y) + mean(dat0$y))
# var(out$y - out$fitted_y)

#
# saemOut$pars[c(1,3,5,7,2,4)]
# my_fit$fit


ggplot(data=out) +
  geom_line(aes(x=x, y=y, col=id)) +
  geom_line(aes(x=x, y=fitted_y, group=id), linetype=2) +
  facet_wrap(~id)

#
# ggplot(data=out) +
#   geom_line(aes(x=x, y=warped_x, col=id))
# ggplot(saemOutDf) +
#   geom_line(aes(x=x, y=saemOutDf$fitted_warped_x, col=id))
#

