# rm(list=ls())
# gc()

library(plyr)
library(dplyr)
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
                                  f_n_knots = 9),
      trace = FALSE)
})

x0 <- seq(0, 1, length=1001)
B <- bs(x0, df = nrow(my_fit$pars$alpha), intercept = TRUE)
matplot(B %*% my_fit$pars$alpha, type="l")

pred_clust <- laply(my_fit$curves, function(crv){
  which.max(crv$sapprox_cluster_membership)})
init_clust <- dat0$clust[dat0$x==0]
table(pred_clust, init_clust)
matplot(t(my_fit$pars_track$sampled_m_track), type="l")


crv <- my_fit$curves[[1]]
matplot(t((crv$scaled_post_amp_mu[1,] * crv$y_scaling_factor +
  t(crv$fitted_base_shapes) * crv$scaled_post_amp_mu[2,])), type="l")


dat0_aug <- ldply(my_fit$curves, function(crv){
  fitted_crv <- as.data.frame(t((crv$scaled_post_amp_mu[1,] * crv$y_scaling_factor +
       t(crv$fitted_base_shapes) * crv$scaled_post_amp_mu[2,])))
  names(fitted_crv) <- paste0("Fitted", 1:ncol(fitted_crv))
  bind_cols(
    data_frame(id = c(crv$curve_id) + 1,
               x = c(crv$x),
               y = c(crv$y),
               pred_clust = which.max(crv$sapprox_cluster_membership)),
    fitted_crv)
  })
dat0_aug$fitted <-
  dat0_aug$Fitted1 * (dat0_aug$pred_clust==1) +
  dat0_aug$Fitted2 * (dat0_aug$pred_clust==2) +
  dat0_aug$Fitted3 * (dat0_aug$pred_clust==3)
dat0_aug <- join(dat0_aug, dat0[,c("x", "y", "clust")])
dat0_aug %>%
  dplyr::filter(id < 50) %>%
  ggplot() +
  geom_line(aes(x=x, y=y, group=id, col=as.factor(clust)),
            show.legend = FALSE)
dat0_aug %>%
  dplyr::filter(id < 101) %>%
  ggplot() +
  geom_line(aes(x=x, y=y, group=id, col=as.character(clust)),
            show.legend = FALSE) +
  geom_line(aes(x=x, y=fitted, group=id, col=as.character(pred_clust)),
            show.legend = FALSE) +
  facet_wrap(~id)
dat0_aug %>%
  dplyr::filter(id < 101) %>%
  ggplot() +
  geom_line(aes(x=x, y=y, group=id, col=as.factor(clust)),
            show.legend = FALSE) +
  geom_line(aes(x=x, y=Fitted1, group=id), alpha=0.5, col="red") +
  geom_line(aes(x=x, y=Fitted2, group=id), alpha=0.5, col="green") +
  geom_line(aes(x=x, y=Fitted3, group=id), alpha=0.5, col="blue") +
  facet_wrap(~id)





