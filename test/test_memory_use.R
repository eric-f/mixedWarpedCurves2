rm(list=ls())
gc()

library(mixedWarpedCurves2)
library(mixedWarpedCurves)
library(ggplot2)
library(tidyverse)
library(plyr)

# dat0 <- readRDS("../data/test_data.rds")
# dat0 <- readRDS("data/test_data.rds")
# dat0 <- readRDS("~/org/lib/mixedWarpedCurves/data/test_data.rds")

# dat1 <- readRDS("~/org/project/Simulations/data/w1-20x1000/data_w1-20x1000-set_001.rds")
# dat2 <- readRDS("~/org/project/Simulations/data/w1-20x1000/data_w1-20x1000-set_002.rds")
# dat3 <- readRDS("~/org/project/Simulations/data/w1-20x1000/data_w1-20x1000-set_003.rds")
# dat1$id <- as.integer(dat1$id) + 100
# dat2$id <- as.integer(dat2$id) + 200
# dat3$id <- as.integer(dat3$id) + 300
# dat0 <- rbind(dat1, dat2, dat3)
# dat0$id <- as.factor(dat0$id)

dp <- readRDS("~/org/project/SESfda-simulation-and-analysis/data/00-03_seal-2010-18_interp_y.rds")
time <- as.numeric(rownames(dp))
dat0 <- as.data.frame(dp) %>%
  rownames_to_column("Time") %>%
  gather(key=id, value=y, -Time) %>%
  mutate(
    x = as.numeric(Time),
    id = as.factor(id))

system.time({
# lineprof({
  set.seed(1)
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
    track_pars = FALSE)
})
#
#
# saemOutDf <- predict_fsim(saemOut)
# saemOutDf$id <- factor(saemOutDf$id, levels = 1:20)
# ggplot(saemOutDf) +
#   geom_line(aes(x=x, y=y, col=id)) +
#   geom_line(aes(x=x, y=fitted_y, group=id), linetype=2) +
#   facet_wrap(~id)

# bctrl <- attr(dat1, "basis.control")

system.time(my_fit <-
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
                                                                   f_n_knots = 5+2)))

# saveRDS(my_fit, file = "my_saem_fit.rds")
plot(my_fit$aux$proposal_sigma_history, type="l")
plot(my_fit$aux$mh_accept_rate_history)
plot(my_fit$pars_track$sigma2_track[-(1:80)], type="l")
plot(my_fit$pars_track$tau_clusters_track[1, -(1:80)], type="l")
plot(my_fit$pars_track$tau_clusters_track[2, -(1:80)], type="l")
plot(my_fit$pars_track$tau_clusters_track[3, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[1, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[2, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[3, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[4, -(1:80)], type="l")
plot(my_fit$pars_track$alpha_track[5, -(1:80)], type="l")
matplot(t(my_fit$pars_track$big_sigma_track)[-(1:80),1], type="l")
matplot(t(my_fit$pars_track$big_sigma_track)[-(1:80),2], type="l")
matplot(t(my_fit$pars_track$big_sigma_track)[-(1:80),4], type="l")
cluster_pred <- laply(my_fit$curves, function(x){x$sapprox_cluster_pred})
heatmap(cluster_pred)

out <- ldply(my_fit$curves, function(x){
  data.frame(id=x$curve_id,
             x=x$x,
             y=x$y,
             warped_x=x$warped_x,
             fitted_y=x$fitted_y,
             cluster=which.max(x$sapprox_cluster_pred))
})
out$id <- as.factor(out$id+1)
out$cluster <- as.factor(out$cluster)
#
# var(dat0$y - dat0$y.shape)
# var(saemOutDf$y * sd(dat0$y) + mean(dat0$y) - saemOutDf$fitted_y * sd(dat0$y) + mean(dat0$y))
# var(out$y - out$fitted_y)

#
# saemOut$pars[c(1,3,5,7,2,4)]
my_fit$pars


ggplot(data=out) +
  geom_line(aes(x=x, y=y, col=id, group=id), show.legend = FALSE) +
  scale_y_reverse() +
  geom_line(aes(x=x, y=fitted_y, group=id), linetype=2, show.legend = FALSE) +
  facet_wrap(~id)

ggplot(data=out) +
  geom_line(aes(x=x, y=y, col=cluster, group=id)) +
  scale_y_reverse() +
  facet_wrap(~cluster)

ggplot(data=out) +
  geom_line(aes(x=x, y=warped_x, group=id, col=cluster))


# ggplot(data=out) +
#   geom_line(aes(x=x, y=warped_x, col=id))
# ggplot(saemOutDf) +
#   geom_line(aes(x=x, y=saemOutDf$fitted_warped_x, col=id))
#

