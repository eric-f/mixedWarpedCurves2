rm(list=ls())
gc()

library(mixedWarpedCurves2)
# library(mixedWarpedCurves)
library(ggplot2)
library(tidyverse)
library(plyr)

dp <- readRDS("~/org/project/SESfda-simulation-and-analysis/data/00-03_seal-2010-18_interp_y.rds")
time <- as.numeric(rownames(dp))
dat0 <- as.data.frame(dp) %>%
  rownames_to_column("Time") %>%
  gather(key=id, value=y, -Time) %>%
  mutate(
    x = as.numeric(Time),
    id = as.factor(id))

# dat0 <- readRDS("data/test_data.rds")

system.time(my_fit <-
              mixedWarpedCurves2::fsim(y = dat0$y,
                                       obs_time = dat0$x,
                                       curve_id = dat0$id,
                                       n_clust = 3,
                                       saem_control = control_saem(n_saem_iter = 100,
                                                                   n_saem_burn = 50,
                                                                   saem_step_seq_pow = 1,
                                                                   n_mcmc_burn = 5,
                                                                   n_core = 1,
                                                                   accept_rate_window = 5,
                                                                   prop_sigma = 1e-2,
                                                                   need_centering = FALSE,
                                                                   accept_rate_lb = 0.17,
                                                                   accept_rate_ub = 0.33,
                                                                   h_n_knots = 1+2,
                                                                   f_n_knots = 1+2)))
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
matplot(t(apply(my_fit$pars_track$kappa_clusters_track[,1,-(1:80)], 2, cumsum)), type="l")
matplot(t(apply(my_fit$pars_track$kappa_clusters_track[,2,-(1:80)], 2, cumsum)), type="l")
matplot(t(apply(my_fit$pars_track$kappa_clusters_track[,3,-(1:80)], 2, cumsum)), type="l")
matplot(t(my_fit$pars_track$sampled_m_track), type="l")


cluster_pred <- laply(my_fit$curves, function(x){x$sapprox_cluster_membership})
heatmap(cluster_pred)


out <- ldply(my_fit$curves, function(x){
  data.frame(id=x$curve_id,
             x=x$x,
             y=x$y,
             warped_x=x$warped_x,
             fitted_y=x$fitted_y,
             cluster=which.max(x$sapprox_cluster_membership),
             pc1=x$sapprox_cluster_membership[1],
             pc2=x$sapprox_cluster_membership[2],
             pc3=x$sapprox_cluster_membership[3])
})
out$id <- as.factor(out$id+1)
out$cluster <- as.factor(out$cluster)

my_fit$pars$tau_clusters

out %>%
  ggplot(data=.) +
  geom_line(aes(x=x, y=y, col=cluster, group=id), show.legend = FALSE) +
  scale_y_reverse() +
  facet_grid(~cluster)

ggplot(data=out) +
  geom_line(aes(x=x, y=warped_x, group=id, col=cluster))



