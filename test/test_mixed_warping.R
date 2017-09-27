rm(list=ls())
gc()

source("R/sim_warping_mixture.R")
library(mixedWarpedCurves2)
library(plyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)

# dp <- readRDS("~/org/project/SESfda-simulation-and-analysis/data/00-03_seal-2010-18_interp_y.rds")
# time <- as.numeric(rownames(dp))
# dat0 <- as.data.frame(dp) %>%
#   rownames_to_column("Time") %>%
#   gather(key=id, value=y, -Time) %>%
#   mutate(
#     x = as.numeric(Time),
#     id = as.factor(id))

# dat0 <- readRDS("data/test_data.rds")
# dat0 <- readRDS("~/org/project/Simulations/data/w1-20x100/data_w1-20x100-set_000.rds")


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
                            tau = 10,
                            mu_sh = -25, mu_sc = 500,
                            sd_sh = 10, sd_sc=50, sd_err = 10)


ggplot(dat0) +
  geom_line(aes(x=x, y=warped_x, group=id, col=clust),
            show.legend = FALSE)

ggplot(dat0) +
  geom_line(aes(x=x, y=y, group=id, col=clust),
            show.legend = FALSE)


system.time(my_fit <-
              mixedWarpedCurves2::fsim_unimodal(
                y = dat0$y,
                obs_time = dat0$x,
                curve_id = dat0$id,
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



mixed_warping_fitted <- ldply(my_fit$curves, function(crv){
  data.frame(id = crv$curve_id,
             x=crv$x, y=crv$y,
             pred_warped_x=crv$warped_x,
             fitted_y=crv$fitted_y,
             cluster=which.max(crv$sapprox_cluster_membership))
})
mixed_warping_fitted$id <- as.factor(mixed_warping_fitted$id)
mixed_warping_fitted$cluster <- as.factor(mixed_warping_fitted$cluster)
laply(my_fit$curves, function(crv){
  crv$sapprox_cluster_membership
})
## Warping coefficients
mixed_warping_coef <- laply(my_fit$curves, function(crv){
  crv$sapprox_w
})
mixed_warping_fitted$cluster_km <-
  as.factor(rep(kmeans(mixed_warping_coef, 3, nstart = 20)$cluster,
                each=201))


# ggplot(dat0) +
#   geom_line(aes(x=x, y=y, group=id, col=clust),
#             show.legend = FALSE)
#
# mixed_warping_fitted %>%
#   ggplot() +
#   geom_line(aes(x=x, y=y, group=id, col=cluster),
#             show.legend = FALSE)
#
#
#
# ggplot(dat0) +
#   geom_line(aes(x=x, y=warped_x, group=id, col=clust),
#             show.legend = FALSE)
#
mixed_warping_fitted %>%
  ggplot() +
  geom_line(aes(x=x, y=pred_warped_x, group=id, col=cluster),
            show.legend = FALSE)
#
#
#
# head(dat0)
# head(mixed_warping_fitted)

fit_all <- cbind.data.frame(dat0, mixed_warping_fitted[,c("pred_warped_x", "fitted_y", "cluster", "cluster_km")])

fit_all %>%
  dplyr::filter(x==0) %>%
  xtabs(~clust+cluster, data=.)

matplot(my_fit$pars_track$sampled_m_track, type="l")
my_fit$curves[[1]]$init_clust

ldply(my_fit$curves, function(crv){
  c(crv$init_clust, which.max(crv$sapprox_cluster_membership))}) %>%
  table

idx <- 200:1000
plot(my_fit$pars_track$sigma2_track[idx], type="l")
plot(my_fit$pars_track$mu_a_track[1,idx], type="l")
plot(my_fit$pars_track$mu_a_track[2,idx], type="l")
matplot(t(my_fit$pars_track$kappa_clusters_track[,1,idx]), type="l")
matplot(t(my_fit$pars_track$kappa_clusters_track[,2,idx]), type="l")
matplot(t(my_fit$pars_track$kappa_clusters_track[,3,idx]), type="l")
