# rm(list=ls())
# gc()
#
library(mixedWarpedCurves2)
# library(plyr)
# library(dplyr)
# library(ggplot2)
# library(tibble)
# library(tidyr)
# library(gtools)
# #
# #
# kappa0 <- c(1,2,2,1)
# kappa0 <- kappa0 / sum(kappa0)
# kappa1 <- c(0,1,1,0)
# kappa1 <- kappa1 / sum(kappa1)
# kappa2 <- c(1.2,1,1,1.2)
# kappa2 <- kappa2 / sum(kappa2)
# kappa3 <- c(0,1,2,1)
# kappa3 <- kappa3 / sum(kappa3)
# kappa4 <- c(1,1,0,0)
# kappa4 <- kappa4 / sum(kappa4)
# #
# dat0 <- sim_warping_mixture(200, rep(1/3, 3),
#                             rbind(kappa1,
#                                   kappa2,
#                                   kappa3),
#                             ni = 1000,
#                             tau = 20,
#                             mu_sh = -25, mu_sc = 500,
#                             sd_sh = 10, sd_sc=50, sd_err = 10)
#
# ggplot(dat0) +
#   geom_line(aes(x=x, y=y, group=id, col=clust), alpha=0.2)
# #
# saveRDS(object = dat0, file = "data/test_mixed_warped_curves.rds")

dat0 <- readRDS("~/org/lib/mixedWarpedCurves2/data/test_mixed_warped_curves.rds")

init_clust <- as.integer(dat0$clust[dat0$x==0]) - 1
system.time(my_fit <-
              try(mixedWarpedCurves2::fsim_mixed_warped_curves(
                y = dat0$y,
                obs_time = dat0$x,
                curve_id = dat0$id,
                init_clust = init_clust,
                n_clust = 3,
                saem_control = control_saem(n_saem_iter = 2000,
                                            n_saem_burn = 100,
                                            saem_step_seq_pow = 1,
                                            n_mcmc_burn = 5,
                                            n_core = 1,
                                            accept_rate_window = 5,
                                            prop_sigma = 1e-2,
                                            need_centering = TRUE,
                                            accept_rate_lb = 0.17,
                                            accept_rate_ub = 0.33,
                                            h_n_knots = 1+4),
                trace = TRUE)))

saveRDS(my_fit, "~/org/lib/mixedWarpedCurves2/test/test_mixed_warped_curves.rds")
quit("no")


rm(list=ls())
library(ggplot2)
library(mclust)
library(fda)
library(dplyr)

dat0 <- readRDS("~/org/lib/mixedWarpedCurves2/data/test_mixed_warped_curves.rds")
init_clust <- as.integer(dat0$clust[dat0$x==0]) - 1
my_fit <- readRDS("~/org/lib/mixedWarpedCurves2/test/test_mixed_warped_curves.rds")


str(my_fit$curves[[1]])
out <- plyr::ldply(my_fit$curves, function(crv){
  data.frame(id = as.character(crv$curve_id+1),
             x = c(crv$x),
             y = c(crv$y),
             warped_x = c(crv$warped_x),
             fitted_y = c(crv$fitted_y),
             pred_clust = LETTERS[which.max(crv$sapprox_cluster_membership)])
})
out$true_clust <- dat0$clust
out$pred_clust <- factor(out$pred_clust, levels = c("A", "B", "C"))

ggplot(out) +
  geom_line(aes(x=x, y=y, group=id),
            show.legend = FALSE) +
  geom_line(aes(x=x, y=fitted_y, color=id),
            show.legend = FALSE) +
  facet_grid(pred_clust~true_clust)

out %>%
  filter(id==as.character(sample(200, 1))) %>%
  ggplot() +
  geom_line(aes(x=x, y=y), col="tomato") +
  geom_line(aes(x=x, y=fitted_y))


pred_clust <- plyr::laply(my_fit$curves, function(crv){which.max(crv$sapprox_cluster_membership)})
table(init_clust, pred_clust)
adjustedRandIndex(init_clust, pred_clust)

out$true_warped_x <- dat0$warped_x
ggplot(out) +
  # geom_line(aes(x=true_warped_x, y=warped_x, group=id)) +
  geom_line(aes(x=x, y=warped_x, group=id)) +
  geom_line(aes(x=x, y=true_warped_x, group=id), col="salmon") +
  facet_grid(pred_clust~true_clust)

sapprox_a <- plyr::laply(my_fit$curves, function(crv){crv$sapprox_a})
(mu_a <- colMeans(sapprox_a))
my_fit$aux$f_order
f0 <- create.bspline.basis(c(0, 1), breaks = c(my_fit$aux$f_break_points))
alpha = my_fit$pars$alpha * mu_a[2] + mu_a[1]
curve(-4*x*(1-x) * 500 - 25, from=0, to=1, col="red")
lines(fd(alpha, f0))
# lines(fd(alpha*mu_a[2] + mu_a[1], f0), col="blue")

prob <- plyr::laply(my_fit$curves, function(crv){crv$sapprox_cluster_membership})

my_fit$pars
my_fit$pars$p_clusters
colSums(my_fit$pars$kappa_clusters)


my_fit$curves[[1]]$sapprox_log_dw

